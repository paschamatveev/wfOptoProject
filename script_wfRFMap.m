%githubDir = 'C:\Users\nadia\Documents\GitHub\wf gui';
githubDir= 'C:\GitHub\wf gui';

addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines')))


%%

mn = 'AL_0035';
td = '2024-08-22';
en = 4;
serverRoot = expPath(mn, td, en);
%% exp numbers: at least the camera folder and the stimulus folder
success = matchBlocks2Timeline(mn,td,[1 4], []);

%% receptive fields 

rfExpNum = 4;
block = getBlockFile(mn, td, rfExpNum);
expRoot = expPath(mn, td, rfExpNum);

pdTimes = readNPY(fullfile(expRoot, 'photodiodeFlips.timestamps.npy'));

[stimTimeInds, stimPositions, stimArray, xPos, yPos] = computeSparseNoiseSignals(block);

stimOnTimes = pdTimes(stimTimeInds{1});
stimOffTimes = pdTimes(stimTimeInds{2});
stimPos = stimPositions{1}; 

%% save rf map stimulus info

writeNPY(stimOnTimes, fullfile(expRoot, 'rfOnTimes.npy'))
writeNPY(stimOffTimes, fullfile(expRoot, 'rfOffTimes.npy'))
writeNPY(stimPos, fullfile(expRoot, 'rfPositions.npy'))