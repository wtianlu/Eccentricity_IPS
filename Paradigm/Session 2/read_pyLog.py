#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 19:19:46 2019

@author: Lulu

190501 Check logfiles, after new script with new timing logging

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#dataFolder = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Paradigm/Session 2/edfpyData'
dataFolder = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data/2019.03 Data/190503_test_original_logging/'
logFname = [f for f in os.listdir(dataFolder) if ".log" in f]
print logFname

fi = int(raw_input('Choose file: '))
with open(dataFolder + os.sep + logFname[fi],'r') as f: 
    logdata = f.readlines()   
#logdata = logdata[109:]
start2end = [l for l in logdata if 'Participant' in l]

trialstarttime = np.array([float(l.split()[0]) for l in logdata if 'TRIAL start' in l])
trialredcross1time = np.array([float(l.split()[0]) for l in logdata if 'show_ShapeStim_redcross 0' in l])
trialcuetime = np.array([float(l.split()[0]) for l in logdata if 'show_BufferImageStim_cue' in l])
trialredcross2time = np.array([float(l.split()[0]) for l in logdata if 'show_ShapeStim_redcross 1' in l])
trialstimtime = np.array([float(l.split()[0]) for l in logdata if 'show_STIM' in l])
trialmasktime = np.array([float(l.split()[0]) for l in logdata if 'show_MASK' in l])
trialredcross3time = np.array([float(l.split()[0]) for l in logdata if 'show_ShapeStim_redcross 2' in l])
trialprobetime = np.array([float(l.split()[0]) for l in logdata if 'show_TextStim_probe' in l])
trialreporttime = np.array([float(l.split()[0]) for l in logdata if 'show_ShapeStim_whitecross - REPORTING' in l])
trialendtime = np.array([float(l.split()[0]) for l in logdata if 'TRIAL end' in l])
trial12start = np.array([float(l.split()[0]) for l in logdata if 'show_white cross start - 12t' in l])
trial12end = np.array([float(l.split()[0]) for l in logdata if 'show_white cross end - 12t' in l])
alltrialtimes = np.array([trialstarttime,trialredcross1time,trialcuetime,trialredcross2time,trialstimtime,\
                          trialmasktime,trialredcross3time,trialprobetime,trialreporttime,trialendtime])
PDtrialtiming = pd.DataFrame({'1.Start':trialstarttime,'2.Cue':trialcuetime,'3.Delay1':trialredcross2time,
                              '4.Stim':trialstimtime,'5.Mask':trialmasktime,'6.Delay2':trialredcross3time,
                              '7.Probe':trialprobetime,'8.Report':trialreporttime,'9.End':trialendtime})    

pulsetimes = np.array([float(l.split()[0]) for l in logdata if 'EmulatedKey: 5' in l])
plt.plot(np.diff(pulsetimes))
plt.title('MR pulses')
plt.show()
#print [l for l in logdata if 'launchScan: start' in l]

for i in range(len(alltrialtimes)):
    if i<len(alltrialtimes)-1:
        tmp = (alltrialtimes[i+1]-alltrialtimes[i])*1000
    else:
        tmp = (alltrialtimes[-1]-alltrialtimes[0])*1000
#    plt.plot(tmp,marker='*')
    plt.hist(tmp)
    plt.title("Dur: "+str(round(np.mean(tmp),3))+" +- "+  str(round(np.std(tmp),3)))
    plt.show()  
print '8 blocks in between of 3.2s (theoretically)'
print trial12end-trial12start

# get cues and probes (NB time of creation, not of presentation!)
datatrial = [int(l.split()[-1]) for l in logdata if 'TRIAL start' in l]
datacue = [int(l[-7]) for l in logdata if 'EXP 	show_BufferImageStim_cue' in l]
dataconf = [int(l[-2]) for l in logdata if 'EXP 	show_BufferImageStim_cue' in l] 
datatarget=[[],[],[],[]]
for i in range(4):
    datatarget[i] = [l[-3] for l in logdata if 'EXP 	TextStim_stim'+str(i)+': text =' in l]
datatarget = map(list,map(None,*datatarget))
dataprobe = [l[-3] for l in logdata if 'EXP 	TextStim_probe: text =' in l]
datachangebool = [0 if dataprobe[i] in datatarget[i] else 1 for i in range(len(dataprobe))]
dataprobeloc = [datatarget[i].index(dataprobe[i]) if dataprobe[i] in datatarget[i] else -1 for i in range(len(dataprobe))]

reportintervalind = [[logdata.index(l) for l in logdata if 'show_TextStim_probe' in l],
                      [logdata.index(l) for l in logdata if 'TRIAL end' in l]]
responselines = [logdata.index(l) for l in logdata if 'Keypress: ' in l]
dataresponses,dataresponsesraw,dataRT,datarespind,datatrialvalidity,datarespcorrect=[],[],[],[],[],[]
responsemap = ['r','b','y']
datarespcorrect=[]
for t in range(len(reportintervalind[0])):
    r = [l[-2] for l in logdata[reportintervalind[0][t]:reportintervalind[1][t]+1] if 'Keypress: ' in l]
    rt = [float(l.split()[0]) for l in logdata[reportintervalind[0][t]:reportintervalind[1][t]+1] if 'Keypress: ' in l]
    dataresponsesraw.append(r)
    dataRT.append(rt-trialprobetime[t])
    datarespind.append([logdata.index(l) for l in logdata[reportintervalind[0][t]:reportintervalind[1][t]+1] if 'Keypress: ' in l])
    # get validity
    if dataprobeloc[t]==-1:v=2
    elif datacue[t]<4 and dataprobeloc[t]==datacue[t]: v=0
    elif datacue[t]==4 and dataprobeloc[t] in [0,3]:v=0
    elif datacue[t]==5 and dataprobeloc[t] in [1,2]:v=0
    else: v=1
    datatrialvalidity.append(v)
    # get correct (1) or not (0), or DK (2) or invalid multiple (3) or invalid no response (4)
    if len(r)==0:
        datarespcorrect.append(4)
        dataresponses.append(-1)
    elif len(r)>1:
        datarespcorrect.append(3)
        dataresponses.append(3)
    else:
        dataresponses.append(responsemap.index(r[0]))
        if responsemap.index(r[0])==2:datarespcorrect.append(2)
        elif responsemap.index(r[0])==datachangebool[t]: datarespcorrect.append(1)
        else: datarespcorrect.append(0)
    
dataNresp = [len(s) for s in dataresponsesraw]    

print 'Valid responses: '+str(sum(dataNresp))+' of total logged keypresses: '+str(len(responselines))

PDtrialdata = pd.DataFrame({'a_Trialind':datatrial,'b_Cue':datacue,'c_Conf':dataconf,'d_Probeloc':dataprobeloc,
                            'e_Probepresent':datachangebool,'f_Validity':datatrialvalidity,'g_targets':datatarget})

sd=7
plannedtrialonsets = range(sd*10,sd*10+32*12,32)+range(sd*10+32*13,sd*10+32*25,32)+\
range(sd*10+32*26,sd*10+32*38,32)+range(sd*10+32*39,sd*10+32*51,32)+range(sd*10+32*52,sd*10+32*64,32)+\
range(sd*10+32*65,sd*10+32*77,32)+range(sd*10+32*78,sd*10+32*90,32)+range(sd*10+32*91,sd*10+32*103,32)
plannedtrialonsets = np.array([float(t)/10 for t in plannedtrialonsets])

plt.plot(trialcuetime-plannedtrialonsets-0.2,marker='o')
    







