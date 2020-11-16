# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 6 - 190103 after merging f and B together, all changes in evernote implemented
190104 using while wait for mask, delay 2, probe letter (doesn't work as well when gaze data is being saved?)
190105 split off from functional: removed references to fMRI and those settings
190109 Added iView calibration and check for existing file names
190204 changed setup to variable x/y-coordinates, fixed radius, constrained theta
190205 Add DK rate to feedback after every block
190206 Moved code shared between session 1 and 2 to separate script
@author: Lena reference Draft6_windows

VERSION 190305 random masks and final extra block
"""

from psychopy import visual, event, core
from random import shuffle
import os,sys
import datetime
import TVAEcc1B_func as tf

if 0:
    from iViewXAPI import *
    res = iViewXAPI.iV_Connect(c_char_p('127.0.0.1'), c_int(4444), c_char_p('127.0.0.1'), c_int(5555))
    if res == 1:
        res = iViewXAPI.iV_GetSystemInfo(byref(systemData))
        print "iV_GetSystemInfo: " + str(res)
        print "Samplerate: " + str(systemData.samplerate)
        print "iViewX Version: " + str(systemData.iV_MajorVersion) + "." + str(systemData.iV_MinorVersion) + "." + str(systemData.iV_Buildnumber)
        print "iViewX API Version: " + str(systemData.API_MajorVersion) + "." + str(systemData.API_MinorVersion) + "." + str(systemData.API_Buildnumber)
        eye=True
        
        @WINFUNCTYPE(None, CSample)
        def sample_callback(sample):
            fsample.write('{timestamp} {left_gazeX} {left_gazeY}\n'.format(
                timestamp=sample.timestamp, 
                left_gazeX=sample.leftEye.gazeX, 
                left_gazeY=sample.leftEye.gazeY))
        
             
        @WINFUNCTYPE(None, CEvent)  
        def event_callback(event):
            fevent.write('{fixstart} {eye} {duration} {posX} {posY}\n'.format(
                fixstart=event.startTime,
                eye=event.eye,
                duration=event.duration, 
                posX=event.positionX,
                posY=event.positionY))  
    else:
        print "Connection failed"
        eye=False
else:
    eye=False

def run_once(f):
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            wrapper.has_run = True
            return f(*args, **kwargs)
    wrapper.has_run = False
    return wrapper

@run_once 
def prepfortrial(t):
    cue_ind = cueseq[t]
    prb_ind = prbseq[t]
    con_ind = conseq[t]
    timing = frames[timseq[t]]
    # select cue
    cue_buffer = cue_list[con_ind][cue_ind]

    # set target letters
    shuffle(letters)
    stiml = letters[0:4]
    if t in B3TRIALS: 
        [stim_display[i].setPos((xcoo[i],ycoo[xyconfigs[conseq[t]][i]])) for i in range(len(xcoo))]
        [mask_display[i].setPos((xcoo[i],ycoo[xyconfigs[conseq[t]][i]])) for i in range(len(xcoo))] # random masks
    [stim_display[i].setText(stiml[i]) for i in range(len(stiml))]
    stim_buffer = visual.BufferImageStim(win,stim=stim_display + [fixation])
    # set masks
#    mask_buffer = mask_display[con_ind][t%len(mask_display[0])]
    shuffle(mask_imgs)
    [mask_display[i].setImage(mask_imgs[i]) for i in range(len(xcoo))]
    mask_buffer = visual.BufferImageStim(win,stim=mask_display + [fixation])
    # find probe letter
    probe_letter = letters[prb_ind]                                 
    probe_display.setText(probe_letter)

    sts = str(t+1)+"\t"+str(cue_ind)+"\t"+ \
    str(prb_ind)+"\t"+str(timseq[t])+"\t"+str(timing)+"\t"+str(con_ind)+"\tT:"+''.join(stiml)+"\tP:"+probe_letter
    
    return cue_buffer,stim_buffer,mask_buffer,probe_display,timing,sts

def pressedesc():
    win.close()
    core.quit()
    if eye:
        fsample.close()
        fevent.close()
        iViewXAPI.iV_Disconnect()
    
################################################################################            
# Set path, create subject file + enter subj details
cwd = os.getcwd()
dataFolder = os.getcwd() + os.path.sep+'smipyData'+os.path.sep
if not os.path.exists(dataFolder): os.makedirs(dataFolder)

try: 
    trialtype = int(raw_input("Enter fPractice [0] or fTask [1]: "))
except: 
    print "Error - Enter 0 or 1"
    sys.exit()
    
if trialtype == 0:
    trial= "BPractice"
    NRUNS,NBLOCKS=1,8
    eye=False
elif trialtype == 1: #Task
    trial= "BTask"
    NRUNS,NBLOCKS = 13,8

subj_id = raw_input("Enter participant ID [xx]: ")
subj_id, fname = tf.loadFnames(subj_id,trial,dataFolder)

# PROCEDURES SETTINGS
NTRIALS,BREAKTRIALS,B3TRIALS,waittimes,letters=tf.getprocsettings(NRUNS,NBLOCKS)
factor,responsekeys, stimsize, masksize, crosssize, E1, E2,xcoo, ycoo, xyconfigs = tf.getdisplaysettings()

# LOAD TRIAL SEQUENCE #
if trialtype==0: 
    cueseq,prbseq,conseq,timseq = tf.loadPseqs()
else:
    cueseq,prbseq,conseq,timseq = tf.loadEseqs(cwd,subj_id)
    
###############################################################################################################################################
# Calibrate eye tracker
fsample,fevent="",""
calibstring = ""
if eye:
    si = raw_input("Calibrate? y/[n]: ")
    while si == 'y':
        iViewXAPI.iV_ContinueEyetracking()
        iViewXAPI.iV_SetupCalibration(byref(calibrationData)) 
        res = iViewXAPI.iV_Calibrate()
        print res
        iViewXAPI.iV_Validate()
        iViewXAPI.iV_GetAccuracy(accuracyData,0)
        calibstring = calibstring + "Start calibration: "+str(datetime.datetime.now())+ "\n" +"Accuracy X: " +\
        str((accuracyData.deviationLX+accuracyData.deviationRX)/2) + " deg\n" +    "Accuracy Y: " + \
        str((accuracyData.deviationLY+accuracyData.deviationRY)/2) + " deg\n"
        
        print "Accuracy X: " + str((accuracyData.deviationLX+accuracyData.deviationRX)/2) + " deg"
        print "Accuracy Y: " + str((accuracyData.deviationLY+accuracyData.deviationRY)/2) + " deg"
        
        si = raw_input("Re-calibrate? y/[n]: ")
    print " "
    print "Finished calibration, switch to DVI"
    calibstring = calibstring + "Finished calibration: "+str(datetime.datetime.now())+ "\n"

################################################################################
# PREPARE VISUAL STIMULI #

#Set up window
win = visual.Window(size=[1920, 1080],color=[-1,-1,-1],monitor='testMonitor',
                             units ='pix', screen=0, winType= "pyglet", fullscr=True)
win.winHandle.activate()
m = event.Mouse(win=win)
m.setVisible(False)

# TIMING SETTINGS
FR = win.getActualFrameRate()
if FR<80:
    frames= [1, 2, 4, 6, 9, 12, 9, 18] #behavioural 0-5, fMRI 6, 300ms 7
    print "!!!WARNING!!! Low window refresh rate detected: " + str(round(FR,2))
else:
    frames= [1, 2, 5, 8, 14, 20, 15, 30]
    
# pre-load cues, masks, stimulus, probe display, ...
#instructions,fbtxt,infoscrtxt,fixation,fixation2,cue_list,mask_display,stim_display,probe_display = tf.loadvisuals(\
#                                                                                win,cwd,masksize,stimsize,crosssize,factor,xcoo,ycoo,xyconfigs,trial)
instructions,fbtxt,infoscrtxt,fixation,fixation2,cue_list,mask_imgs,mask_display,stim_display,probe_display = tf.loadvisuals(\
                                                                                win,cwd,masksize,stimsize,crosssize,factor,xcoo,ycoo,xyconfigs,trial)

###############################################################################
# Write to file
file=open(dataFolder+fname, "w") #create participant file (rn it's one file per phase)
file.writelines(["Participant "+subj_id, '\t', trial,'\n'])
file.write("Start: "+str(datetime.datetime.now())+"\nFrame rate: " + str(round(FR)) + "\nVisual stimuli (E5 pixels): " + str(stimsize))
file.write("\n\n" + calibstring)
file.write("\n\n#\tCueI\tPrbI\tTimI\tFrames\tConfig\tTargets\tProbe\tOnset\tActTime\tResp\tRT\tTrialtimings\n")
file.close()

if eye:    
    fsample = open(dataFolder+fname[0:-4]+'sample.txt','w')
    fevent  = open(dataFolder+fname[0:-4]+'event.txt','w') 
    res = iViewXAPI.iV_SetSampleCallback(sample_callback)
    res = iViewXAPI.iV_SetEventCallback(event_callback)
    iViewXAPI.iV_PauseEyetracking()

#################################################################################################################################
#################################################################################################################################
# START EXPERIMENT
#################################################################################################################################
#################################################################################################################################

current_trial = tf.infoscreen(win,instructions,infoscrtxt,E1,E2,xcoo,ycoo,crosssize,responsekeys,eye,fsample,fevent)

# define clocks
gnclock = core.Clock()
runclock = core.Clock()
# pre-prepare first trial displays
cue_buffer,stim_buffer,mask_buffer,probe_display,timing,sts = prepfortrial(current_trial)
prepfortrial.has_run = False

while current_trial<NTRIALS:
    ################################################################################
    # PREPARING THIS TRIAL'S DISPLAYS
    gnclock.reset() 
    trialtiming,responsetiming = [],[]
    initial_string = "-" #reset                     
    
    ###########################################################################
    # check progress through the trials: break blocks, sub blocks...
    write = False
    if current_trial in BREAKTRIALS:
        correct,incorrect,dks=0,0,0
        fbtxt.setText("Start block "+str(current_trial/96+1)+" of " + str(NRUNS))
        fbtxt.draw() 
        win.flip()

        waittime = waittimes[0]
        write,wait = True,True
        while wait:
             for key in event.getKeys():
                 if key in ['5','return','space']:
                    wait=False
                    fixation2.draw()  
                    win.flip()
                    if eye: iViewXAPI.iV_ContinueEyetracking()
                    event.clearEvents()
                    runclock.reset() 
                    gnclock.reset()   
                 elif key in ['escape']:
                     pressedesc()
        
    elif current_trial in B3TRIALS: 
        waittime = waittimes[1] #every 12 trials
    else:
        waittime=0.2
    
    if write and eye: savetxtstring = str(datetime.datetime.now()) +"\n"+sts
    elif write: savetxtstring = str(datetime.datetime.now()) +"\n"+sts
    else: savetxtstring = sts

    ################################################################################
    # SHOW THIS TRIAL'S DISPLAYS    
    crossflip = True
    while 1:
        if gnclock.getTime()>waittime-1/FR:        
            break
        elif waittime-gnclock.getTime() < 0.2:    
            fixation.draw()
            win.flip()
        else:
            fixation2.draw()
            win.flip()
            
    if eye: fsample.write('\nTrial {ct} onset {ts}\n'.format(ct=current_trial+1, ts = runclock.getTime()))
    
    savetxtstring = savetxtstring+"\t"+str(round(runclock.getTime(),5))
  
    trialtiming.append(gnclock.getTime()*1000) #2 start cue
    # cue of 300ms
    for i in range(frames[7]):
        cue_buffer.draw()
        win.flip()

    trialtiming.append(gnclock.getTime()*1000) #3 start del
    # second delay of 300ms
    for i in range(frames[7]):
        fixation.draw()
        win.flip()

    #display stimuli for x ms
    trialtiming.append(gnclock.getTime()*1000) #4    start stim
    for i in range(timing):
        stim_buffer.draw()
        win.flip()
    
    trialtiming.append(gnclock.getTime()*1000) #5 start mask
    
    if eye: fsample.write('Mask {ct} onset {ts}\n'.format(ct=current_trial+1, ts = runclock.getTime()))
    
    #masks for 300ms
    for i in range(frames[7]):
        mask_buffer.draw()
        win.flip()
    
    trialtiming.append(gnclock.getTime()*1000) #6 start delay
    
    # third delay of 200ms
    for i in range(frames[5]):
        fixation.draw()
        win.flip()

    trialtiming.append(gnclock.getTime()*1000) #7 start probe
    #probe letter of 200ms
    for i in range(frames[5]):
        probe_display.draw()
        win.flip()
    trialtiming.append(gnclock.getTime()*1000) #8 start response
    
    fixation2.draw() #accounted for in response time
    win.flip()
        
    ################################################################################
    # GET RESPONSES
    if eye: fsample.write('Report {ct} onset {ts}\n'.format(ct=current_trial+1, ts = runclock.getTime()))

    while  1:
        if not prepfortrial.has_run and current_trial<NTRIALS-1: 
            cue_buffer,stim_buffer,mask_buffer,probe_display,timing,sts = prepfortrial(current_trial+1)
        if gnclock.getTime()*1000-trialtiming[-1]>=1500: 
            break
        for key,kt in event.getKeys(timeStamped=gnclock):
            if key in ['escape']:
                pressedesc()
            elif key in responsekeys:
                responsetiming.append(kt*1000)
                initial_string = initial_string + str(responsekeys.index(key)%3+1)           
            else:
                initial_string = initial_string+'('+key+')'
                responsetiming.append(kt*1000)

    prepfortrial.has_run = False
           
    # Update performance
    prb_ind = prbseq[current_trial]
    if len(initial_string)>1:
        if (prb_ind < 4 and initial_string[1]=='1') or (prb_ind>3 and initial_string[1]=='2'):
            correct+=1
        elif (prb_ind < 4 and initial_string[1]=='2') or (prb_ind>3 and initial_string[1]=='1'):
            incorrect+=1
        elif initial_string[1]=='3':
            dks += 1

    # Write trial timing information
    trialtiming.append(gnclock.getTime()*1000) #9 end response
    savetxtstring = savetxtstring + "\t" + str(round(trialtiming[3]-trialtiming[2],1))+"\tR"+initial_string+"\t"+\
    str([round(elem,1) for elem in responsetiming])+"\tTT:"+ str([round(elem, 1) for elem in trialtiming]) +"\n"

    #show feedback
    if current_trial+1 in BREAKTRIALS:
        if eye: iViewXAPI.iV_PauseEyetracking()
        if (correct+incorrect)==0: incorrect+=1 #to avoid ZeroDivisionError
        showtxt = "End of block "+str(current_trial/96+1)+"\n\nPercentage correct: " + str(correct*100/(correct+incorrect)) + "%"
        if current_trial == NTRIALS-1:
            showtxt = showtxt + "\n\n\nEnd of the " + trial[1:].lower() +trialtype*" \n\n Thank you for your participation!"
        elif current_trial == NTRIALS-96:
            showtxt = showtxt + "\n\n\n Next: unguided run\nPlease call the experimenter" 
        fbtxt.setText(showtxt)
        fbtxt.draw() 
        win.flip()
        savetxtstring = savetxtstring + showtxt+"\nPercentage don't know: " + str(round(dks/.96)) + '%\n'+str(datetime.datetime.now()) +"\n"
        
        wait = True
        while wait:
             for key in event.getKeys():
                 if key in ['5','return','space']:
                    wait=False
                 elif key in ['escape']:
                     file=open(dataFolder+fname, "a+")
                     file.write(savetxtstring)
                     file.close()   
                     pressedesc()
                 elif key in ['m','d']:
                     win.fullscr = False           
                     win.winHandle.set_fullscreen(False)
                     win.winHandle.minimize()
                     win.flip()
                    
                     cont = raw_input('Press enter to continue with the next block')
                     win.winHandle.maximize()
                     win.winHandle.activate()
                     win.fullscr=True
                     win.setColor([-1,-1,-1])
                     win.winHandle.set_fullscreen(True)
                     fbtxt.setText("Press enter to proceed")
                     fbtxt.draw() 
                     win.flip()
    
    file=open(dataFolder+fname, "a+")
    file.write(savetxtstring)
    file.close()   
    current_trial +=1


pressedesc()




