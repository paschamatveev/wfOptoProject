%% Script to analyze a new widefield+opto recording

githubDir = ['C:\Users\nadia\Documents\GitHb\wf gui'];
% githubDir= 'C:\GitHub\wf gui';

%% get galvos, laserpowers, etc as before - assuming step function, and putting all in #1 folder

addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir , 'Pipelines'))) % steinmetzlab/Pipelines
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
% addpath(genpath(fullfile(githubDir, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis

mn = 'AL_0033'; 
td = '2024-12-17';
ca_en = 1; % widefield

serverRoot = expPath(mn, td, ca_en);

expHz = [false false]; %


%% process - the original way, not multiple exps, and WITHOUT 40hz
% check timeline signals

gx = readNPY(fullfile(serverRoot,'galvoX.raw.npy'));
gy = readNPY(fullfile(serverRoot,'galvoY.raw.npy'));
laser = readNPY(fullfile(serverRoot,'lightCommand.raw.npy'));
ts = readNPY(fullfile(serverRoot,'lightCommand.timestamps_Timeline.npy'));

%
t = tsToT(ts, numel(laser));

[~, laserOn, laserOff] = schmittTimes(t, laser, [0.025 0.075]);
laserPower = round(interp1(t, laser, laserOn+0.005), 1); % interp halfway through stim when laser has been on
galvoXPos = round(interp1(t, gx, laserOn), 1);
galvoYPos = round(interp1(t, gy, laserOn), 1);

galvoPos = [galvoXPos, galvoYPos];

% face camera timestamps
cameraTrigger = readNPY(fullfile(serverRoot, 'cameraTrigger.raw.npy'));
cameraTriggerTL = readNPY(fullfile(serverRoot, 'cameraTrigger.timestamps_Timeline.npy'));

t =  interp1(cameraTriggerTL(:, 1), cameraTriggerTL(:, 2), 1:numel(cameraTrigger))';
[flipTimes, flipsUp, flipsDown] = schmittTimes(t, cameraTrigger, [1 4]);

expTimes = readNPY(fullfile(serverRoot, 'expStartStop.timestamps_Timeline.npy'));
expStartStop = readNPY(fullfile(serverRoot, 'expStartStop.raw.npy'));
[times,expStart,expEnd] = schmittTimes(t,expStartStop,[.5 .5]); %i think these are in seconds?
sampPerSec= ts(2)/ts(2,2); % sample/set

%% write

% write laser on, off, power, galvo
writeNPY(laserPower, fullfile(serverRoot, 'laserPowers.npy'));
writeNPY(galvoXPos, fullfile(serverRoot, 'galvoXPositions.npy'));
writeNPY(galvoYPos, fullfile(serverRoot, 'galvoYPositions.npy'));
writeNPY(laserOn, fullfile(serverRoot,'laserOnTimes.npy'));
writeNPY(laserOff,fullfile(serverRoot,'laserOffTimes.npy'));

[positions, ~, posLabels] = unique(galvoPos, 'rows');
[powers, ~, powerLabels] = unique(laserPower);
 
% write cam
writeNPY(flipsUp, fullfile(serverRoot, 'cameraFrameTimes.npy'));

%% process - saves all to one folder, and depends on you hard-coding what experiments are 40hz and which are not

%read in raw signals 
gx = readNPY(fullfile(serverRoot,'galvoX.raw.npy'));
gy = readNPY(fullfile(serverRoot,'galvoY.raw.npy'));
laser = readNPY(fullfile(serverRoot,'lightCommand.raw.npy'));
ts = readNPY(fullfile(serverRoot,'lightCommand.timestamps_Timeline.npy'));
expTimes = readNPY(fullfile(serverRoot, 'expStartStop.timestamps_Timeline.npy')); 
expStartStop = readNPY(fullfile(serverRoot, 'expStartStop.raw.npy')); 

t = tsToT(ts, numel(laser)); %seconds timing
sampPerSec= ts(2)/ts(2,2); % sample/sec 
[times,expStart,expEnd] = schmittTimes(t,expStartStop,[2 4.5]); %find times in s of exp start/stop

% finding laser ons and offs times for a 40hz signal
sigName = 'lightCommand';
[tt, v] = getTLanalog(mn, td, ca_en, sigName);

tInd = 1;
traces(tInd).t = tt; %tsToT(ts, numel(raw), explicit times
traces(tInd).v = v; %raw times
traces(tInd).name = sigName;

threshold=.15;
stimlen=0.05;
stimTimes = tt(v(2:end)>threshold & v(1:end-1)<=threshold); %find times where light is on, has to be at least 2V
stimOffsets = tt(v(2:end)<threshold & v(1:end-1)>=threshold);

ds = find(diff([0;stimTimes])>stimlen); % find indeces where there is a stim, time between has to be greater than .05
stimStarts = stimTimes(ds); %ds

gapDur = [stimTimes; max(tt)]-[0;stimOffsets];
stimEnds = stimOffsets(gapDur(2:end)>stimlen);

%find the powers for each stim
pwRnd = zeros(length(stimEnds),1); 
for i = 1:length(stimEnds)
    startpt=stimStarts(i);
    endpt=stimEnds(i);
    
    %find duration
    dur = round(endpt-startpt,2);
    if dur == 0.02 
        dur = 0.025;
    end

    if dur < 0.025 %if it finds a duration that is less than the minimum (0.025), it skips
        continue
    else
        samples = (dur) * sampPerSec; %end sec - start sec * (samples/sec)
        interp=zeros(samples,1); 
        interp(:,1) = interp1(tt,laser,linspace(startpt,endpt,samples)); %find the laser values 
        pwsMax = max(interp(:,1)); % find max of the laser                  between the start and end sec
        pwRnd(i) = round(pwsMax/2,1);

        % for j = 1:length(expStart)
        %     if stimEnds(i) > expStart(j) && stimEnds(i) < expEnd(j) %if its a 40hz exp, divide by 2. if not, 
        %         if expHz(j) == true                                 % its a step, so the max is the right value
        %             pwRnd(i) = round(pwsMax/2,1);
        %         else 
        %             pwRnd(i) = round(pwsMax,1);
        %         end
        %     end
        % end
    end
end

galvoXPos = round(interp1(t, gx, stimStarts), 1); 
galvoYPos = round(interp1(t, gy, stimStarts), 1);

galvoPos = [galvoXPos, galvoYPos];

% face camera timestamps
cameraTrigger = readNPY(fullfile(serverRoot, 'cameraTrigger.raw.npy'));
cameraTriggerTL = readNPY(fullfile(serverRoot, 'cameraTrigger.timestamps_Timeline.npy'));

t =  interp1(cameraTriggerTL(:, 1), cameraTriggerTL(:, 2), 1:numel(cameraTrigger))';
[flipTimes, flipsUp, flipsDown] = schmittTimes(t, cameraTrigger, [1 4]);


%% write

% write laser on, off, power, galvo
writeNPY(stimStarts, fullfile(serverRoot, 'laserOnTimes.npy'));
writeNPY(stimEnds, fullfile(serverRoot, 'laserOffTimes.npy'));
writeNPY(pwRnd, fullfile(serverRoot, 'laserPowers.npy'));
writeNPY(galvoXPos, fullfile(serverRoot, 'galvoXPositions.npy'));
writeNPY(galvoYPos, fullfile(serverRoot, 'galvoYPositions.npy'));

[positions, ~, posLabels] = unique(galvoPos, 'rows');
[powers, ~, powerLabels] = unique(pwRnd);
 
% write cam
writeNPY(flipsUp, fullfile(serverRoot, 'cameraFrameTimes.npy'));


%% separating exps + laser powers + saving into different folders
% this one is quite broken
expTimes = readNPY(fullfile(serverRoot, 'expStartStop.timestamps_Timeline.npy'));
expStartStop = readNPY(fullfile(serverRoot, 'expStartStop.raw.npy'));
[times,expStart,expEnd] = schmittTimes(t,expStartStop,[2 4.5]); %i think these are in seconds?
sampPerSec= ts(2)/ts(2,2); % sample/set
%prep to find laser powers
sigName = 'lightCommand';
[tt, v] = getTLanalog(mn, td, ca_en, sigName);

tInd = 1;
traces(tInd).t = tt; %tsToT(ts, numel(raw), explicit times
traces(tInd).v = v; %raw times
traces(tInd).name = sigName;

for i = 1:size(expStart)
    serverRoot_exp = expPath(mn, td, i+1);

    %convert to samples 
    indStart = find(t==expStart(i));
    indEnd = find(t==expEnd(i));
    
    %find lasers, galvos, times, camera,etc during the experiment
    laser_exp = laser(indStart:indEnd); %expStart and expStop are in SECONDS. i need the INDEX where that SECOND happens. 
    t_exp = t(indStart:indEnd);
    gx_exp = gx(indStart:indEnd);
    gy_exp = gy(indStart:indEnd);
    cameraTrigger_exp = cameraTrigger(indStart:indEnd);
    tt_exp = tt(indStart:indEnd);
    v_exp = v(indStart:indEnd);
    
    %find the laser times 
    threshold=0.05;
    stimlen=0.05;
    stimTimes = tt_exp(v_exp(2:end)>threshold & v_exp(1:end-1)<=threshold); %find times where light is on, has to be at least 2V
    stimOffsets = tt_exp(v_exp(2:end)<threshold & v_exp(1:end-1)>=threshold);
    
    ds = find(diff([0;stimTimes])>stimlen); % find indeces where there is a stim, time between (i think) has to be greater than .05
    stimStarts = stimTimes(ds); %ds
    
    gapDur = [stimTimes; max(tt)]-[0;stimOffsets];
    stimEnds = stimOffsets(gapDur(2:end)>stimlen);

    pwsMax=zeros(length(stimEnds),1); 
    pwsRnd = zeros(length(stimEnds),1); 
    
    %find laser powers
    for j = 1:length(stimEnds)
        startpt=stimStarts(j);
        endpt=stimEnds(j);
    
        dur = round(endpt-startpt,2);
        if dur == 0.02 % putting a bandaid on a boat leak right here
            dur = 0.025;
        end
        samples = (dur) * sampPerSec; %end sec - start sec * (samples/sec)

        interp=zeros(samples,1);
        interp(:,1) = interp1(t_exp,laser_exp,linspace(startpt,endpt,samples));
        pwsMax(j) = max(interp(:,1));

        % check if the experiment is oscillating or not
        % if it is oscillating - halve the max. if not, don't.
        %hard-coding in the right values because rounding is weird
        if expHz(i) == true
            pw = round(pwsMax(j)/2,1);
            pwsRnd(j) = pw;
        else
            pw = round(pwsMax(j),1);
            pwsRnd(j) = pw;
        end
    end
    %get galvos, flips
    galvoXPosexp = round(interp1(t_exp, gx_exp, stimStarts), 1);
    galvoYPosexp = round(interp1(t_exp, gy_exp, stimStarts), 1);
    
    galvoPos_exp = [galvoXPosexp, galvoYPosexp];
    
    [flipTimes, flipsUp_exp, flipsDown] = schmittTimes(t_exp, cameraTrigger_exp, [1 4]);
    
    writeNPY(stimStarts, fullfile(serverRoot_exp, 'laserOnTimes.npy'));
    writeNPY(stimEnds, fullfile(serverRoot_exp, 'laserOffTimes.npy'));
    writeNPY(pwsRnd, fullfile(serverRoot_exp, 'laserPowers.npy'));
    writeNPY(galvoXPosexp, fullfile(serverRoot_exp, 'galvoXPositions.npy'));
    writeNPY(galvoYPosexp, fullfile(serverRoot_exp, 'galvoYPositions.npy'));
    writeNPY(flipsUp_exp, fullfile(serverRoot_exp, 'cameraFrameTimes.npy')); %might be broken
end

%% TESTING GROUND
t = tsToT(ts, numel(laser));
expTimes = readNPY(fullfile(serverRoot, 'expStartStop.timestamps_Timeline.npy'));
expStartStop = readNPY(fullfile(serverRoot, 'expStartStop.raw.npy'));
[times,expStart,expEnd] = schmittTimes(t,expStartStop,[2 4.5]); %i think these are in seconds?

%find samples per sec
sampPerSec= ts(2)/ts(2,2); % sample/set

% test one exp at a time
trialsec=[2,4]; %setting up how many sec per trial per exp
trEnd = [1,0]; %premaking an array where to put the number of the last trial for each exp
expStartInd = zeros(1, length(expStart));
expEndInd = zeros(1, length(expEnd));
for i = 1:length(expStart)
    sec = expEnd(i)-expStart(i); %find the number of seconds in the exp
    trs = sec/trialsec(i); %find how many trial in the exp
    trEnd(i+1) = round(trs - 1) + trEnd(i); %find number of the last trial in the exp across all
end                                            % the trials that exist 

writeNPY(trEnd,fullfile(serverRoot,'expTrials.npy'))
%% write the appropriate stimstarts, ends, pws, galvos for each exp
% we find activity from the stimstarts/ends, which are in seconds, but are
% the length of the number of trials
% so can use the trial numbers to find the appropriate starts, ends, pws 
for i = 1:length(expStart)
    serverRoot_exp = expPath(mn, td, i+1);
    startTr=trEnd(i);
    endTr = trEnd(i)+ trEnd(i+1);
    disp(i)
    disp(startTr)
    disp(endTr)
    disp(serverRoot_exp)
    % writeNPY(stimStarts(startTr:endTr), fullfile(serverRoot_exp, 'laserOnTimes.npy'));
    % writeNPY(stimEnds(startTr:endTr), fullfile(serverRoot_exp, 'laserOffTimes.npy'));
    % writeNPY(pwsRnd(startTr:endTr), fullfile(serverRoot_exp, 'laserPowers.npy'));
    % writeNPY(galvoXPosexp(startTr:endTr), fullfile(serverRoot_exp, 'galvoXPositions.npy'));
    % writeNPY(galvoYPosexp(startTr:endTr), fullfile(serverRoot_exp, 'galvoYPositions.npy'));

end
%%
%isolate stimends and stimstarts for this exp
sigName = 'lightCommand';
[tt, v] = getTLanalog(mn, td, ca_en, sigName);

tInd = 1;
traces(tInd).t = tt; %tsToT(ts, numel(raw), explicit times
traces(tInd).v = v; %raw times
traces(tInd).name = sigName;

tt_exp = tt(indStart:indEnd);
v_exp = v(indStart:indEnd);
laser_exp=laser(indStart:indEnd);

threshold=0.05;
stimlen=0.05;
stimTimes = tt_exp(v_exp(2:end)>threshold & v_exp(1:end-1)<=threshold); %find times where light is on, has to be at least 2V
stimOffsets = tt_exp(v_exp(2:end)<threshold & v_exp(1:end-1)>=threshold);

ds = find(diff([0;stimTimes])>stimlen); % find indeces where there is a stim, time between (i think) has to be greater than .05
stimStarts = stimTimes(ds); %ds

gapDur = [stimTimes; max(tt)]-[0;stimOffsets];
stimEnds = stimOffsets(gapDur(2:end)>stimlen);

%storage for checking
pwsMean=zeros(length(stimEnds),1); 
pwsHalved=zeros(length(stimEnds),1);
samplesList=zeros(length(stimEnds),1);
durList=zeros(length(stimEnds),1);
lengthInterp = zeros(length(stimEnds),1);
pwsMax=zeros(length(stimEnds),1); 
pwsUnhalved=zeros(length(stimEnds),1); 
pwsRnd = zeros(length(stimEnds),1); 


% find the powers
for i = 1:length(stimEnds)
    startpt=stimStarts(i);
    endpt=stimEnds(i);

    dur = round(endpt-startpt,2);
    if dur == 0.02 % putting a bandaid on a boat leak right here 
        dur = 0.025;
    end
    samples = (dur) * sampPerSec; %end sec - start sec * (samples/sec)
    samplesList(i)=samples;
    durList(i) = dur;

    interp=zeros(samples,1);
    interp(:,1) = interp1(tt_exp,laser_exp,linspace(startpt,endpt,samples));
    pwsMax(i) = max(interp(:,1));

    % check if the experiment is oscillating or not
    % if it is oscillating - halve the max. if not, don't.
    if expHz(exp) == true
        pw = round(pwsMax(i)/2,1);
    else
        pw = round(pwsMax(i),1);
    end
end

disp("done")

%%
h(1)=plot(tt,v,'cyan');
hold on
h(2)=scatter(stimTimes,ones(size(stimTimes)) + 1,'red');
hold on
h(3)=scatter(stimStarts,ones(size(stimStarts)),12,'blue',"*");
hold on
h(4)=scatter(stimEnds,ones(size(stimEnds)),'green');
hold off
legend(h,"raw stimulus (v)",'stimTimes with threshold', 'stimStarts via ds', 'stimEnds via ds')
title("stim threshold " + threshold + " stim len " + stimlen + ' original code')
% xlim([200,250])
% ylim([2.8,3])
saveas(gcf,'stimtimes15.jpg')

%% linearity test

linTestExpNum = 2;

block = getBlockFile(mn,td,linTestExpNum);
success = matchBlocks2Timeline(mn,td, [1 2],[]); %[main exp folder, visresp folder]

stimContrast = [block.paramsValues.Contrast];
stimDur = [block.paramsValues.Duration];
% stimOnTimes = block.events.stim_onTimes;

expRoot = expPath(mn, td, linTestExpNum);

pdTimes = readNPY(fullfile(expRoot, 'photodiodeFlips.timestamps.npy'));
stimOnTimes = pdTimes(1:2:end);
stimOffTimes = pdTimes(2:2:end); 

% save stimuli info
writeNPY(stimContrast, fullfile(expRoot, 'linTestContrasts.npy'));
writeNPY(stimDur, fullfile(expRoot, 'linTestDurations.npy'))
writeNPY(stimOnTimes, fullfile(expRoot, 'linTestOnTimes.npy'));
writeNPY(stimOffTimes, fullfile(expRoot, 'linTestOffTimes.npy'));
