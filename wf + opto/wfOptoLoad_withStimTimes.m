%% Script to analyze a new widefield+opto recording

githubDir = 'C:\Users\nadia\Documents\GitHub\wf gui';
% githubDir= 'C:\GitHub\wf gui';

%% get galvos, laserpowers, etc as before - assuming step function, and putting all in #1 folder

addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir , 'Pipelines'))) % steinmetzlab/Pipelines
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
% addpath(genpath(fullfile(githubDir, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis

mn = 'AB_0032'; 
td = '2024-07-26';
ca_en = 1; % widefield

serverRoot = expPath(mn, td, ca_en);

expHz = [false,false];
%% process - not multiple exps, and without 40hz
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

%% write

% write laser on, off, power, galvo
writeNPY(laserOn, fullfile(serverRoot, 'laserOnTimes.npy'));
writeNPY(laserOff, fullfile(serverRoot, 'laserOffTimes.npy'));
writeNPY(laserPower, fullfile(serverRoot, 'laserPowers.npy'));
writeNPY(galvoXPos, fullfile(serverRoot, 'galvoXPositions.npy'));
writeNPY(galvoYPos, fullfile(serverRoot, 'galvoYPositions.npy'));

[positions, ~, posLabels] = unique(galvoPos, 'rows');
[powers, ~, powerLabels] = unique(laserPower);

% write cam
writeNPY(flipsUp, fullfile(serverRoot, 'cameraFrameTimes.npy'));

%% find exps

expTimes = readNPY(fullfile(serverRoot, 'expStartStop.timestamps_Timeline.npy'));
expStartStop = readNPY(fullfile(serverRoot, 'expStartStop.raw.npy'));
[times,expStart,expEnd] = schmittTimes(t,expStartStop,[2 4.5]); %i think these are in seconds?

%% save galvos, powers, etc for multiple exps
% CURRENTLY ONLY WORKS FOR ONE DURATION. will integrate multiple durations
% once i figure that out!

for i = 1:size(expStart)
    serverRoot_exp = expPath(mn, td, i+1);
    
    %convert to samples 
    indStart = find(t==expStart(i));
    indEnd = find(t==expEnd(i));
    
    laser_exp = laser(indStart:indEnd); %expStart and expStop are in SECONDS. i need the INDEX where that SECOND happens. 
    t_exp = t(indStart:indEnd);
    gx_exp = gx(indStart:indEnd);
    gy_exp = gy(indStart:indEnd);
    cameraTrigger_exp = cameraTrigger(indStart:indEnd);
    
    [~, laserOnexp, laserOffexp] = schmittTimes(t_exp, laser_exp, [0.025 0.075]);
    laserPowerexp = round(interp1(t_exp, laser_exp, laserOnexp+0.005), 1); % interp halfway through stim when laser has been on
    galvoXPosexp = round(interp1(t_exp, gx_exp, laserOnexp), 1);
    galvoYPosexp = round(interp1(t_exp, gy_exp, laserOnexp), 1);
    
    galvoPos_exp = [galvoXPosexp, galvoYPosexp];

    [flipTimes, flipsUp, flipsDown] = schmittTimes(t, cameraTrigger_exp, [1 4]);

    writeNPY(laserOnexp, fullfile(serverRoot_exp, 'laserOnTimes.npy'));
    writeNPY(laserOffexp, fullfile(serverRoot_exp, 'laserOffTimes.npy'));
    writeNPY(laserPowerexp, fullfile(serverRoot_exp, 'laserPowers.npy'));
    writeNPY(galvoXPosexp, fullfile(serverRoot_exp, 'galvoXPositions.npy'));
    writeNPY(galvoYPosexp, fullfile(serverRoot_exp, 'galvoYPositions.npy'));
    writeNPY(flipsUp, fullfile(serverRoot_exp, 'cameraFrameTimes.npy'));
end

%% finding laser ons and offs for a 40hz signal

%nick code
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

ds = find(diff([0;stimTimes])>stimlen); % find indeces where there is a stim, time between (i think) has to be greater than .05
stimStarts = stimTimes(ds); %ds

gapDur = [stimTimes; max(tt)]-[0;stimOffsets];
stimEnds = stimOffsets(gapDur(2:end)>stimlen);

%% separating exps + laser powers + saving into different folders

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
        if dur == 0.02 % putting a bandaid on a boat leak right here (even tho it didnt do anything)
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

%% testing ground to find pows and lengths
%older code from above
t = tsToT(ts, numel(laser));

%find samples per sec
sampPerSec= ts(2)/ts(2,2); % sample/set

% test one exp at a time

exp = 2;
indStart = find(t==expStart(exp));
indEnd = find(t==expEnd(exp));

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