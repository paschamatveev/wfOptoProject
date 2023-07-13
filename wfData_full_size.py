import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import scipy
class wfData:
    '''
    creates all necessary items for working with widefield imaging data
    '''
    def __init__(self, pathSubject, pathStim):
        '''
        creates all the items
        '''
        serverPath = Path(pathSubject)
        self.timeFile = serverPath / 'cameraFrameTimes.npy'
        self.frameTimes = np.squeeze(np.load(self.timeFile))[::2]
        self.svdTemp = np.load(serverPath / 'corr/svdTemporalComponents_corr.npy')
        self.svdSpat = np.load(serverPath / 'blue/svdSpatialComponents.npy')
        self.meanImage = np.load(serverPath / 'blue/meanImage.npy')
        self.stimOnTimes = np.squeeze(np.load(pathStim / 'linTestOnTimes.npy'))
        self.stimContrasts = np.squeeze(np.load(pathStim / 'linTestContrasts.npy'))
        self.stimDurations = np.squeeze(np.load(pathStim / 'linTestDurations.npy'))
        self.getTempComp = scipy.interpolate.interp1d(self.frameTimes, self.svdTemp, axis=0)
        self.spatial = (self.svdSpat).reshape(560*560, -1)
        self.contrasts = np.unique(self.stimContrasts)
        self.durations = np.unique(self.stimDurations)
    def getContrasts():
        return self.contrasts
    def getDurations():
        return self.durations
    def trial_timeFunc(self, start, end, step, oneTime=False, loc=-1):
        self.oneTime = oneTime
        if oneTime==True:
            startTime = self.stimOnTimes[loc] - start
            endTime = self.stimOnTimes[loc] + end
            self.trial_time = np.linspace(startTime, endTime, step)
            self.trial_time = np.array(self.trial_time)
        else:
            self.trial_time = [np.linspace(i-0.1, i+0.4, 100) for i in self.stimOnTimes]
            self.trial_time = np.array(self.trial_time)
        return self.trial_time
    def getTempCompFunc(self, trial_time_input):
        '''
        uses the interpolation function to return values
        '''
        return self.getTempComp(trial_time_input)
    def trial_activityFunc(self, trial_time_input):
        self.trial_activity = self.getTempCompFunc(trial_time_input)
        self.trial_activity = np.array(self.trial_activity)
        return self.trial_activity
    def createVideo(self, trial_activity):
        if self.oneTime==False:
            trial_activity_re = np.mean(trial_activity, axis=0)
        else:
            trial_activity_re = trial_activity
        spatial = (self.svdSpat).reshape(560*560, -1)
        self.video = spatial @ trial_activity_re.T
        self.video = (self.video).reshape(560,560,-1)
        return self.video
    def plotBrain(self, rows=1, cols=1, avg=False, figsizeCol=2, figsizeRow=2):
        if avg==True:
            f = plt.figure(figsize=(cols*figsizeCol, rows*figsizeRow))
            avg_trial_activity = np.mean(self.trial_activity, axis=1)
            self.video = np.mean(self.video, axis=2)
            plt.imshow(self.video[:, :], clim = np.percentile(self.video, (2, 99.9)))
            plt.colorbar()
        else:
            f = plt.figure(figsize=(cols*figsizeCol, rows*figsizeRow))
            gs = mpl.gridspec.GridSpec(rows, cols)

            for i in range(rows*cols):
                ax = plt.subplot(gs[i])
                plt.imshow(self.video[:, :, i*2], clim = np.percentile(self.video, (2, 99.9)))
                plt.colorbar()

            f.tight_layout()
    def createContrastVideos(self, start, stop, step):
        videos=[]
        for contrast in self.contrasts:
            theseIndexes = np.squeeze(np.argwhere(self.stimContrasts == contrast))

            stim = self.stimOnTimes[theseIndexes]
            self.trial_time = [np.linspace(i-start, i+stop, step) for i in stim]

            trial_activity = self.getTempComp(self.trial_time)

            reshape = np.mean(trial_activity, axis=0)
            video = self.spatial @ reshape.T
            video = video.reshape(560, 560, -1)
            videos.append(video)

        videos = np.array(videos)
        return videos
    def createDurationVideos(self, start, stop, step):
        videos=[]
        for dur in self.durations:
            theseIndexes = np.squeeze(np.argwhere(self.stimDurations==dur))
            stim = self.stimOnTimes[theseIndexes]
            self.trial_time=[np.linspace(i-start,i+stop,step) for i in stim]
            trial_activity=self.getTempComp(self.trial_time)
            reshape = np.mean(trial_activity,axis=0)
            video=self.spatial @ reshape.T
            video = video.reshape(560,560,-1)
            videos.append(video)
        videos=np.array(videos)
        return videos
    def createConDurVid(self, start, stop, step, durCount):
        ''' 
        broken
        '''
        videos = [[] for x in range(durCount)]
        for count,duration in enumerate(self.durations[:durCount]):
            theseIndexesDur = np.squeeze(np.argwhere(self.stimDurations == duration))
            for contrast in self.contrasts:
                theseIndexesCon = np.squeeze(np.argwhere(self.stimContrasts == contrast))
                indexes = np.intersect1d(theseIndexesCon, theseIndexesDur)
                stim = self.stimOnTimes[indexes]
                self.trial_time = [np.linspace(i-start, i+stop, step) for i in stim]

                trial_activity = self.getTempComp(self.trial_time)

                reshape = np.mean(trial_activity, axis=0)
                video = self.spatial @ reshape.T
                video = video.reshape(560, 560, -1)
                videos[count].append(video)
                
        videos = np.array(videos)
        return videos
    def corrections(self, frames):
        videoCorr = np.mean(frames, axis=(0, 1))
        minActivity = min(videoCorr)
        if minActivity < 0:
            videoCorr += minActivity*-1
        f0 = np.mean(videoCorr[:6])
        videoCorr = [(x-f0)/(f0+1) for x in videoCorr]
        return np.array(videoCorr)