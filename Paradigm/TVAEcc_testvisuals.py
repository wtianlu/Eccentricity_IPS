# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 6 - 190103 after merging f and B together, all changes in evernote implemented
190104 using while wait for mask, delay 2, probe letter (doesn't work as well when gaze data is being saved
190105 split off from functional: removed references to fMRI and those settings
190109 Added iView calibration and check for existing file names
@author: Lena reference Draft6_windows

"""

from psychopy import visual, event, core
from random import shuffle
import os
import TVAEccfunc_190221 as tf

cwd = os.getcwd()


################################################################################
# PREPARE VISUAL STIMULI #
fMRI = True
if fMRI:
    trial="f"  
    session = 2
else:
    trial="B"
    session = 1

NRUNS = 8
NBLOCKS = 8

"""
dimpx[0] = stimulus size
dimpx[1] = mask size
dimpx[2] = fixation cross size
dimpx[3] = dot size
dimpx[4] = E1
dimpx[5] = E2
"""

################################################################################
#Set up window
r=0.75
win = visual.Window(size=[1920*r, 1080*r],color=[-1,-1,-1],monitor='testMonitor',
                             units ='pix', screen=0, winType= "pyglet", fullscr=False)

NTRIALS,BREAKTRIALS,FBTRIALS,B14TRIALS,B3TRIALS,trialtiming,correct,incorrect,dks,initial_string,letters=tf.getprocsettings(NRUNS,NBLOCKS)
waittimes,factor,responsekeys, stimsize, masksize, crosssize, E1, E2,xcoo, ycoo, xyconfigs = tf.getdisplaysettings(session,0.2)
instructions,fbtxt,infoscrtxt,fixation,fixation2,cue_list,mask_imgs,mask_display,stim_display,probe_display = tf.loadvisuals(win,cwd,masksize,stimsize,crosssize,factor,xcoo,ycoo,trial)

shuffle(letters) #random Ts in every trial
stiml = letters[0:4]
[stim_display[i].setText(stiml[i]) for i in range(len(stiml))]

shuffle (mask_imgs)
[mask_display[i].setImage(mask_imgs[i]) for i in range(len(xcoo))] 
                                   
probe_letter = letters[0]    
probe_display.setText(probe_letter)


# %%
################################################################################


wait = True

infoscrtxt = visual.TextStim(win=win,units='norm', height = 0.1, wrapWidth=1.7,
    alignHoriz='center', text=  "[D]ots\n[C]ue\n[F]ixation\n[S]tim\n[M]ask\n[P]robe\n[B]ack\n[escape]\n[Z] screenshot" + trial[1:].lower())
infoscrtxt.draw()
win.flip()
img=0
whc = 0
wc = -1
wf = 0
fixationlist = [fixation,fixation2]
while wait:
    for key in event.getKeys():
        if key in ['escape']:
            win.close()
            core.quit()
            wait=False
        elif key in ['d']:
            # validation part
            win.color=[-0.5]*3
            win.flip()
            win.flip()
            calibpos = [[0,0],[E2,0],[0,E1],[-E2,0],[0,-E1],[E2,-E1],[E2,E1],[-E2,E1],[-E2,-E1],[0,0]]
            dot = visual.Circle(win=win, radius=0.15*crosssize, fillColor='white',pos=calibpos[0])
            for i in range( len(calibpos)):
                dot.setPos(calibpos[i])
                dot.draw()
            win.flip()
        elif key in ['b']:
            win.color=[-1]*3
            win.flip()
            infoscrtxt.draw()
            win.flip()
            
        elif key in ['c']:
            wc+=1
            cue = cue_list[whc%len(xyconfigs)][wc%len(cue_list[0])]
            print 'cue '+str(wc%len(cue_list[0]))
            [elem.draw() for elem in cue]
            #cue.draw()
            win.flip()
            
        elif key in ['f']:
            wf+=1
            fixationlist[wf%2].draw()
            win.flip()
        elif key in ['s']:
            [stim_display[i].setPos((xcoo[i],ycoo[xyconfigs[whc%len(xyconfigs)][i]])) for i in range(4)]
            whc+=1
            [elem.draw() for elem in stim_display]
            #fixation.draw()
            win.flip()
            print 'config '+str(whc%len(xyconfigs))
        elif key in ['m']:
            [mask_display[i].setPos(stim_display[i].pos) for i in range(4)]
            [elem.draw() for elem in mask_display]
            fixation.draw()
            win.flip()
        elif key in ['p']:
            probe_display.draw()
            win.flip()
        elif key in ['z']:
            img+=1
            win.getMovieFrame()
            win.saveMovieFrames('screenshot'+str(img)+'.png')
            print 'screenshot'+str(img)+'.png'

win.close()







