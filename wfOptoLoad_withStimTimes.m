%% Script to analyze a new widefield+opto recording

githubDir = 'C:\Users\nadia\Documents\GitHub\wf gui';
% githubDir= 'C:\GitHub\wf gui';

%% paths 

addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir , 'Pipelines'))) % steinmetzlab/Pipelines
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
% addpath(genpath(fullfile(githubDir, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis

mn = 'ZYE_0077'; 
td = '2024-07-22';
ca_en = 1; % widefield

serverRoot = expPath(mn, td, ca_en);

%% check timeline signals

gx = readNPY(fullfile(serverRoot,'galvoX.raw.npy'));
gy = readNPY(fullfile(serverRoot,'galvoY.raw.npy'));
laser = readNPY(fullfile(serverRoot,'lightCommand.raw.npy'));
ts = readNPY(fullfile(serverRoot,'lightCommand.timestamps_Timeline.npy'));

%%
t = tsToT(ts, numel(laser));

[~, laserOn, laserOff] = schmittTimes(t, laser, [0.025 0.075]);
laserPower = round(interp1(t, laser, laserOn+0.005), 1); % interp halfway through stim when laser has been on
galvoXPos = round(interp1(t, gx, laserOn), 1);
galvoYPos = round(interp1(t, gy, laserOn), 1);

galvoPos = [galvoXPos, galvoYPos];

writeNPY(laserOn, fullfile(serverRoot, 'laserOnTimes.npy'));
writeNPY(laserOff, fullfile(serverRoot, 'laserOffTimes.npy'));
writeNPY(laserPower, fullfile(serverRoot, 'laserPowers.npy'));
writeNPY(galvoXPos, fullfile(serverRoot, 'galvoXPositions.npy'));
writeNPY(galvoYPos, fullfile(serverRoot, 'galvoYPositions.npy'));

[positions, ~, posLabels] = unique(galvoPos, 'rows');
[powers, ~, powerLabels] = unique(laserPower);

%% face camera timestamps
cameraTrigger = readNPY(fullfile(serverRoot, 'cameraTrigger.raw.npy'));
cameraTriggerTL = readNPY(fullfile(serverRoot, 'cameraTrigger.timestamps_Timeline.npy'));

t =  interp1(cameraTriggerTL(:, 1), cameraTriggerTL(:, 2), 1:numel(cameraTrigger))';
[flipTimes, flipsUp, flipsDown] = schmittTimes(t, cameraTrigger, [1 4]);

writeNPY(flipsUp, fullfile(serverRoot, 'cameraFrameTimes.npy'));

%% find exps

expTimes = readNPY(fullfile(serverRoot, 'expStartStop.timestamps_Timeline.npy'));
expStartStop = readNPY(fullfile(serverRoot, 'expStartStop.raw.npy'));
[times,expStart,expEnd] = schmittTimes(t,expStartStop,[2 4.5]);

%% save exps
sampPerSec= ts(2)/ts(2,2); % sample/sec

for i = 1:size(expTimes)
    serverRoot_exp = expPath(mn, td, i+1);
    
    %find the starts and stops for this experiment and get the right laser
    start_stamp = expStart(i)*sampPerSec;
    end_stamp= expEnd(i)*sampPerSec;

    laser_exp = laser(start_stamp:end_stamp);
    ts_exp = ts(start_stamp:end_stamp);
    cameraTriggerTL_exp = cameraTriggerTL(start_stamp:end_stamp);
    cameraTrigger_exp = cameraTrigger(start_stamp:end_stamp);
    
    t = tsToT(ts_exp, numel(laser_exp));

    %remake galvos, laser times, camera
    [~, laserOn, laserOff] = schmittTimes(t, laser_exp, [0.025 0.075]);
    laserPower = round(interp1(t, laser_exp, laserOn+0.005), 1); % interp halfway through stim when laser has been on
    galvoXPos = round(interp1(t, gx, laserOn), 1);
    galvoYPos = round(interp1(t, gy, laserOn), 1);
    galvoPos = [galvoXPos, galvoYPos];
    t =  interp1(cameraTriggerTL_exp(:, 1), cameraTriggerTL_exp(:, 2), 1:numel(cameraTrigger_exp))';
    [flipTimes, flipsUp, flipsDown] = schmittTimes(t, cameraTrigger_exp, [1 4]);

    writeNPY(laserOn, fullfile(serverRoot_exp, 'laserOnTimes.npy'));
    writeNPY(laserOff, fullfile(serverRoot_exp, 'laserOffTimes.npy'));
    writeNPY(laserPower, fullfile(serverRoot_exp, 'laserPowers.npy'));
    writeNPY(galvoXPos, fullfile(serverRoot_exp, 'galvoXPositions.npy'));
    writeNPY(galvoYPos, fullfile(serverRoot_exp, 'galvoYPositions.npy'));
    writeNPY(flipsUp, fullfile(serverRoot_exp, 'cameraFrameTimes.npy'));
end

%% finding starts and ends
% find starts and ends

%nick code
sigName = 'lightCommand';
[tt, v] = getTLanalog(mn, td, ca_en, sigName);

tInd = 1;
traces(tInd).t = tt; %tsToT(ts, numel(raw), explicit times
traces(tInd).v = v; %raw times
traces(tInd).name = sigName;

threshold=2;
stimlen=0.05;
stimTimes = tt(v(2:end)>threshold & v(1:end-1)<=threshold); %find times where light is on, has to be at least 2V
stimOffsets = tt(v(2:end)<threshold & v(1:end-1)>=threshold);

ds = find(diff([0;stimTimes])>stimlen); % find indeces where there is a stim, time between (i think) has to be greater than .05
stimStarts = stimTimes(ds); %ds

gapDur = [stimTimes; max(tt)]-[0;stimOffsets];
stimEnds = stimOffsets(gapDur(2:end)>stimlen);

writeNPY(stimStarts, fullfile(serverRoot, 'laserOnTimes_test.npy'));
writeNPY(stimEnds, fullfile(serverRoot, 'laserOffTimes_test.npy'));

% plot
% h(1)=plot(tt,v,'cyan');
% hold on
% h(2)=scatter(stimTimes,ones(size(stimTimes)) + 1,'red');
% hold on
% h(3)=scatter(stimStarts,ones(size(stimStarts)),12,'blue',"*");
% hold on
% h(4)=scatter(stimEnds,ones(size(stimEnds)),'green');
% hold off
% legend(h,"raw stimulus (v)",'stimTimes with threshold', 'stimStarts via ds', 'stimEnds via ds')
% title("stim threshold " + threshold + " stim len " + stimlen + ' original code')
% xlim([200,250])
% ylim([2.8,3])
% saveas(gcf,'stimtimes15.jpg')

%% find pws and lengths
%older code from above
t = tsToT(ts, numel(laser));

%find samples per sec
sampPerSec= ts(2)/ts(2,2); % sample/sec

%one power onset per stimEnd/Start
pwsMean=zeros(length(stimEnds),1); 
pwsRnd=zeros(length(stimEnds),1);

for i = 1:length(stimEnds)
    startpt=stimStarts(i);
    endpt=stimEnds(i);

    samples = (endpt-startpt) * sampPerSec; %end sec - start sec * (samples/sec)
    samples= round(samples);

    interp=zeros(samples,1);
    interp(:,1) = interp1(t,laser,linspace(startpt,endpt,samples)); %over time, take laser values, assign to values from start:end:samples

    pwsMean(i) = mean(interp); %saved for testing
    pwsRnd(i) = round(pwsMean(i),1);
end

writeNPY(pwsRnd, fullfile(serverRoot, 'laserPowers_test.npy'));
