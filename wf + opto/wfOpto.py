import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import scipy

import pytoolsAL as ptAL

class wfOpto:
    '''
    contains some functions for common operations with wf-opto data
    '''
    def __init__(self, pathSubject, listExps=None):
        '''
        pathSubject - path to the experiments
        listExps - a list of lists. used if there are multiple experiments done in one session. 
                    each sublist is the indexes of the corresponding experiment. 
        '''
        serverPath = Path(pathSubject)
        self.timeFile = serverPath / 'cameraFrameTimes.npy'
        self.frameTimes = np.squeeze(np.load(self.timeFile))[::2] # every other frame - we want blue only
        self.svdTemp = np.load(serverPath / 'corr/svdTemporalComponents_corr.npy')
        self.svdSpat = np.load(serverPath / 'blue/svdSpatialComponents.npy')
        self.svdSpatFull = self.svdSpat[:,:,:500]
        
        self.meanImage = np.load(serverPath / 'blue/meanImage.npy')
        self.laserOn = np.squeeze(np.load(serverPath / 'laserOnTimes.npy'))
        self.laserOff = np.squeeze(np.load(serverPath / 'laserOffTimes.npy'))
        self.laserPowers = np.squeeze(np.load(serverPath /'laserPowers.npy'))
        self.galvoX = np.squeeze(np.load(serverPath/'galvoXPositions.npy'))
        self.galvoY = np.squeeze(np.load(serverPath/'galvoYPositions.npy'))
        self.px, self.py, self.ncomps = self.svdSpatFull.shape
        
        self.svdSpat = self.svdSpatFull.reshape(self.px*self.py, self.ncomps)
        self.tToWf = scipy.interpolate.interp1d(self.frameTimes, self.svdTemp, axis=0, fill_value='extrapolate')
        self.spatial = self.svdSpatFull.reshape(560*560,-1)
        if listExps == None:
            listExps = np.array([np.arange(len(self.laserPowers))])
        self.listExps = listExps
        self.pulseLengths = []
        for count,time in enumerate(self.laserOff):
            length = self.laserOff[count]-self.laserOn[count]
            length = round(length,2)
            self.pulseLengths.append(length)
        self.pulseLengths = np.array(self.pulseLengths)
    def tToWFManual(self, time):
        self.trial_activity = self.tToWf(time)
        return self.trial_activity
    def oneTrial(self, start, stop, step, trial=0, exp=0):
        '''
        creates video for a single trial
        user must employ python indexing (starting with 0) 
        '''
        startTime = self.laserOn[self.listExps[exp][trial]] + start
        endTime = self.laserOn[self.listExps[exp][trial]] + stop
        
        trial_time = np.linspace(startTime, endTime, step)
        trial_activity = self.tToWf(trial_time)
        
        dwf = [np.diff(i, prepend=i[0]) for i in trial_activity.T]
        dwf = np.array(dwf)
        
        spatial = self.svdSpatFull.reshape(560*560, -1)
        video = self.spatial @ dwf
        video = video.reshape(560, 560, -1)
        
        n_cols = 5
        n_rows = 10
        f = plt.figure(figsize=(n_cols*2, n_rows*2))
        gs = mpl.gridspec.GridSpec(n_rows, n_cols)
        for i in range(50):
            ax = plt.subplot(gs[i])
            plt.imshow(video[:, :, i*2],clim = np.percentile(video, (2, 99.9)), cmap='bwr')
            plt.title(f't = {trial_time[i*2]:.2f}s')
            plt.colorbar()
            
        f.tight_layout()
    def allTrials(self, start, stop, step,exp=0):
        '''
        videos for all trials
        '''
        trial_time_all = [np.linspace(i+start, i+stop, step) for i in self.laserOn[self.listExps[exp]]]
        trial_activity_all = self.tToWf(trial_time_all)
        trial_activity_all = np.mean(trial_activity_all, axis=0)
        
        dwf = [np.diff(i, prepend=i[0]) for i in trial_activity_all]
        dwf = np.array(dwf)
        
        spatial = (self.svdSpatFull).reshape(560*560, -1)
        video = spatial @ dwf.T
        video = video.reshape(560, 560, -1)
        
        n_cols = 5
        n_rows = 10
        f = plt.figure(figsize=(n_cols*2, n_rows*2))
        gs = mpl.gridspec.GridSpec(n_rows, n_cols)
        for i in range(50):
            ax = plt.subplot(gs[i])
            brain=plt.imshow(video[:, :, i*2], clim = np.percentile(video, (2, 99.9)), cmap='bwr')
            plt.colorbar(brain)
            
        f.tight_layout()
    def fullAvg(self,start,stop,step,trials=None,exp=0):
        '''
        creates one image for average of all trials that are selected using the trials argument 
        '''
        if trials.any() == None:
            trials = np.linspace(self.listExps[exp][0],self.listExps[exp][-1], self.listExps[exp][-1]-self.listExps[exp][0], dtype=int)
        trial_time_all = [np.linspace(i+start, i+stop, step) for i in self.laserOn[self.listExps[exp][trials]]]
        trial_activity_all = self.tToWf(trial_time_all)
        trial_activity_all = np.mean(trial_activity_all, axis=0)
        
        dwf = [np.diff(i, prepend=i[0]) for i in trial_activity_all]
        dwf = np.array(dwf)
        
        avg_trial_activity = np.mean(dwf, axis=1)
        
        spatial = self.svdSpatFull.reshape(560*560, -1)
        videoAvg = spatial @ dwf.T
        videoAvg = videoAvg.reshape(560,560,-1)
        videoAvg = np.mean(videoAvg, axis=2)
        plt.imshow(videoAvg[:,:], clim = np.percentile(videoAvg, (2, 99.9)), cmap='bwr')
    def compareAvgs(self, trials=None, start=0, stop=100, n_col=10, n_row=10, exp=0):
        '''
        creates image of avg activity for all trials between start and stop
        uses trials argument to decide which trials to do.
        user has to know how they want their trials to be represented, and how their trials will fit in their desired rows/cols
        '''
        allVideos = []
        if trials.any()== None:
            trials = np.linspace(0,100,100,dtype=int)
        for trial in trials:
            startTime = self.laserOn[self.listExps[exp][trial]]
            endTime = self.laserOn[self.listExps[exp][trial]] + .5
            
            trial_time = np.linspace(startTime, endTime, 100)
            trial_activity = self.tToWf(trial_time)
            dwf = [np.diff(i, prepend=i[0]) for i in trial_activity.T]
            dwf = np.array(dwf)
            
            avg_trial_activity = np.mean(dwf, axis=1)
        
            videoAvg = self.spatial @ avg_trial_activity.T
            videoAvg = videoAvg.reshape(560,560,-1)
            videoAvg = np.mean(videoAvg, axis=2)
        
            allVideos.append(videoAvg)

        allVideos = np.array(allVideos)

        f = plt.figure(figsize=(n_col*3, n_row*3))
        gs = mpl.gridspec.GridSpec(n_row, n_col)
        clim = np.percentile(allVideos,(2,99.9))
        
        for count,trial in enumerate(trials):
            ax = plt.subplot(gs[count])
            thisVideo = allVideos[count]
            plt.imshow(thisVideo, clim=clim, cmap='bwr')
            
            ax = ptAL.plotting.apply_image_defaults(ax)
            plt.title("trial " + str(trial + 1))
            x = self.galvoX[self.listExps[exp][trial]]
            y = self.galvoY[self.listExps[exp][trial]]
            length = self.pulseLengths[self.listExps[exp][trial]]
            power = self.laserPowers[self.listExps[exp][trial]]
            plt.text(0,700,f'position: [{x},{y}] \n length: {length:.5f} \n power: {power}', fontsize=10)
            plt.xlabel([])
            plt.ylabel([])
            cb = ptAL.plotting.add_colorbar(ax)
            
        f.tight_layout()
    def trackPixel(self, x, y, start=-.2,stop=.7,step=100,exp=0,trialStart=None,trialStop=None):
        '''
        pixel activity over desired trials
        '''
        if trialStop == None:
            trialStart = self.listExps[exp][0]-self.listExps[exp][0]
            trialStop = (self.listExps[exp][-1]-1) - self.listExps[exp][0]
        timeScale = np.linspace(start,stop,step)
        for trial in range(trialStart,trialStop):
            startTime = self.laserOn[self.listExps[exp][trial]] + start
            endTime = self.laserOn[self.listExps[exp][trial]] + stop
            
            trial_time = np.linspace(startTime, endTime, 100)
            trial_activity = self.tToWf(trial_time)
            
            dwf = [np.diff(i, prepend=i[0]) for i in trial_activity.T]
            dwf = np.array(dwf)
            
            spatial = self.svdSpatFull.reshape(560*560, -1)
            video = spatial @ dwf # can multiply by spatial indexing by just one pixel 
            video = video.reshape(560, 560, -1)
                
            plt.plot(timeScale, video[420,450],marker='.',c='mediumorchid')
        plt.xlabel('Trial Time (milisec)')
        plt.ylabel('Activity')
        plt.title(f'Activity of pixel {x}, {y} over trials {trialStart} - {trialStop}')
    def standardError(self, x, y, exp=0, trials=None,start=-.2,stop=.7,step=100,brain=False):
        '''
        standard error for a certain pixel over certain trials
        '''
        if trials.any() == None:
            trials = np.arange(self.listExps[exp][0]-self.listExps[exp][0],(self.listExps[exp][-1]-1) - self.listExps[exp][0], 1)
        xval = len(trials)
        pixelVals = np.zeros((xval,100))
        spatial = self.svdSpatFull.reshape(560*560, -1)
        for trial in trials:
            startTime = self.laserOn[self.listExps[exp][trial]] - .1
            endTime = self.laserOn[self.listExps[exp][trial]] + .7
            
            trial_time = np.linspace(startTime, endTime, 100)
            trial_activity = self.tToWf(trial_time)
        
            dwf = [np.diff(i, prepend=i[0]) for i in trial_activity.T]
            dwf = np.array(dwf)
        
            video = spatial @ trial_activity.T
            video = video.reshape(560, 560, -1)
        
            for timePt in range(100):
                # fills in activity of the desired pixel for all trials. 333x100, is filling in the 100 for the trial out of 333. uses 100 timepts. 
                pixelVals[trial][timePt] = video[x,y,timePt]
            timeScale = np.linspace(start,stop,step)
        if brain:
            trial_time_all = [np.linspace(i+start, i+stop, step) for i in self.laserOn[self.listExps[exp]]]
            trial_activity_all = self.tToWf(trial_time_all)
            trial_activity_all = np.mean(trial_activity_all, axis=0)
            
            dwf = [np.diff(i, prepend=i[0]) for i in trial_activity_all]
            dwf = np.array(dwf)
            
            avg_trial_activity = np.mean(dwf, axis=1)
            
            videoAvg = self.spatial @ dwf.T
            videoAvg = videoAvg.reshape(560,560,-1)
            videoAvg = np.mean(videoAvg, axis=2)

            fig, (ax1,ax2) = plt.subplots(1,2)
            ax1.plot(timeScale,np.mean(pixelVals[:,:], axis=0),color='darkslateblue')
            ax1.fill_between(timeScale, \
                             np.mean(pixelVals[:,:], axis=0)-scipy.stats.sem(pixelVals[:,:],axis=0),\
                             np.mean(pixelVals[:,:],axis=0)+scipy.stats.sem(pixelVals[:,:],axis=0), color='lavender')
            ax1.set_xlabel("time in milisec")
            ax1.set_ylabel("activity")
            ax1.set_title(f'Standard error for pixel {x}, {y}')

            ax2.imshow(videoAvg[:,:], clim = np.percentile(videoAvg, (2, 99.9)), cmap='bwr')
            ax2.scatter([x],[y],color='green',marker='v')
        else:
            plt.plot(timeScale,np.mean(pixelVals[:,:], axis=0),color='darkslateblue')
            plt.fill_between(timeScale, \
                             np.mean(pixelVals[:,:], axis=0)-scipy.stats.sem(pixelVals[:,:],axis=0),\
                             np.mean(pixelVals[:,:],axis=0)+scipy.stats.sem(pixelVals[:,:],axis=0), color='lavender')
            plt.xlabel("time in milisec")
            plt.ylabel("activity")
            plt.title(f'Standard error for pixel {x}, {y}')