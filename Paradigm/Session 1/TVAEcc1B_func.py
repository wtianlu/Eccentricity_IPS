#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 15:20:01 2019

@author: Lulu

TVAeccentricity study functions

getprocsettings(NRUNS,NBLOCKS)
    return NTRIALS,BREAKTRIALS,FBTRIALS,B14TRIALS,B3TRIALS,trialtiming,correct,incorrect,dks, initial_string,letters
getdisplaysettings(session)
    return waittimes,factor,responsekeys, stimsize, masksize, crosssize, E1, E2, xcoo, ycoo, xyconfigs
loadFnames(session,subj_id,trial)

loadPseqs(session)
loadEseqs(cwd,subj_id,session)
loadvisuals(win,cwd,masksize,stimsize,crosssize,factor,E1,E2,trial,pos=4)
loadguidedots(win,E1,E2,crosssize)
infoscreen(session,testbuttons,getcurtrial,wait,win,instructions,infoscrtxt,E1,E2,crosssize,responsekeys,eye=False,fname_eye='')

VERSION 190305 random masks
"""
import os
import math
from psychopy import visual, event, core
from random import shuffle
if 0: from iViewXAPI import *
        
def getprocsettings(NRUNS,NBLOCKS):
    nt = 12 #trials per block
    
    NTRIALS = nt*NBLOCKS*NRUNS
    BREAKTRIALS = range(0,NTRIALS+1,NTRIALS/NRUNS)
    B3TRIALS = range(0,NTRIALS,NTRIALS/(NRUNS*NBLOCKS))

    waittimes = [3.4,3.4]
    letters = ["A", "B", "D", "E", "F", "G", "H", "J", "K", "L","M", "N",
          "O", "P", "R", "S", "T", "V", "X","Z"]
    
    return NTRIALS,BREAKTRIALS,B3TRIALS,waittimes,letters


def getdisplaysettings():
    
    factor = [1.617, 1, 1,1.617,0.3923] # including S and N: 1.65
    sizes = [2.24, 2.77, 0.4, 5, 10]

    dimpx = [round(math.tan(math.radians(d))*600/0.2767) for d in sizes]
    responsekeys = ['1','2','3','num_1','num_2','num_3']
    stimsize,masksize,crosssize,E1,E2 = dimpx

    i = 6
    th_lim = [math.pi/i, math.pi/i]
    xcoo = [dimpx[-1]*math.cos(math.pi-th_lim[1]/2),dimpx[-2]*math.cos(math.pi-th_lim[0]/2),dimpx[-2]*math.cos(th_lim[0]/2),dimpx[-1]*math.cos(th_lim[1]/2)]
    ycoo = [dimpx[-1]*math.sin(math.pi-th_lim[1]/2),dimpx[-2]*math.sin(math.pi+th_lim[0]/2),dimpx[-2]*math.sin(math.pi-th_lim[0]/2),dimpx[-1]*math.sin(math.pi+th_lim[1]/2)]
    xyconfigs = [[0,1,1,0],[3,2,2,3]] #top to bottom
    
    return factor,responsekeys, stimsize, masksize, crosssize, E1, E2, xcoo, ycoo, xyconfigs



################################################################################################################################################################
#
#   LOAD PARTICIPANT SEQUENCES STUFF
#
################################################################################################################################################################

def loadFnames(subj_id,trial,dataFolder):
    if int(subj_id)<10 and len(subj_id)<2:
        subj_id = '0'+subj_id
    fname = "P" + subj_id + "_" + trial + ".txt"
    
    while os.path.isfile(dataFolder + fname):
        print "File "+fname+" already exists"
        si = raw_input("Enter new (full) filename or press enter to overwrite: ")
        if len(si)>1:
            fname = si + fname.endswith(".txt")*".txt"            
        else:
            break
    print "New "+fname+" created"
    return subj_id, fname



def loadPseqs():
    """behavioural practice (taken from random run of p1sequence_B; timing is just shuffled for 2nd half)"""
    cueseq = [5,1,1,4,5,1,4,4,4,5,5,1,3,6,3,2,2,6,2,6,6,3,2,3,5,4,4,1,4,5,1,4,5,1,5,1,3,2,3,3,6,6,2,3,6,2,2,6,              
          1,1,4,5,1,1,4,5,4,5,5,4,3,6,3,3,2,6,3,6,2,2,6,2,5,5,4,5,4,1,5,1,1,1,4,4,2,6,6,6,3,3,2,3,2,2,3,6]
    prbseq = [7,1,5,4,5,7,7,5,4,1,4,1,2,7,3,7,4,5,2,2,3,5,5,7,4,4,6,2,4,6,8,8,8,2,2,6,6,2,3,3,4,6,6,8,2,8,2,8,              
          1,1,4,7,7,5,1,5,7,1,4,5,3,1,3,5,7,3,7,7,5,2,5,2,4,8,8,1,4,6,6,1,1,8,6,4,2,3,6,2,3,3,2,8,8,6,6,8]
#    prbseq = [7,1,5,3,5,7,7,5,4,2,4,1,3,7,2,7,4,5,2,2,3,5,5,7,4,4,6,1,2,6,8,8,8,4,1,6,6,1,3,3,4,6,6,8,2,8,2,8,              1,3,4,7,7,5,1,5,7,1,4,5,3,1,4,5,7,3,7,7,5,2,5,2,3,8,8,1,4,6,6,1,2,8,6,4,3,3,6,2,3,1,2,8,8,6,6,8]
    cueseq = [i-1 for i in cueseq]
    prbseq = [i-1 for i in prbseq]
    conseq = [0]*24+[1]*24+[0]*24+[1]*24
    timseq = [5,5,5,5,5,5,5,5,5,5,5,5,4,5,5,5,5,4,5,5,4,4,4,4,4,5,5,5,5,4,5,5,4,4,4,4,4,5,5,4,3,4,4,3,3,5,5,3,         
              4,5,5,4,3,4,4,3,3,5,5,3,4,5,5,4,3,4,4,3,3,5,5,3,0,3,0,4,5,2,5,4,2,3,1,1,0,3,0,4,5,2,5,4,2,3,1,1]
    return cueseq,prbseq,conseq,timseq

def loadEseqs(cwd,subj_id):
    cueseq,prbseq,conseq,timseq=[],[],[],[]
    f = open(cwd + '/sequences/p' + subj_id + 'sequence_B.txt', 'r')
    tmp = f.readlines()
    f.close()  

    for i in range(len(tmp)):
        cueseq.append(int(tmp[i].split('\t')[0]))
        prbseq.append(int(tmp[i].split('\t')[1]))
        timseq.append(int(tmp[i].split('\t')[2]))
        conseq.append(int(tmp[i].split('\t')[3]))
    if min(cueseq)==1:
        cueseq = [i-1 for i in cueseq]
    if min(prbseq)==1:
        prbseq = [i-1 for i in prbseq]
    if min(conseq)==1:
        conseq = [i-1 for i in conseq]
    if min(timseq)==1:
        timseq = [i-1 for i in timseq]

    return cueseq,prbseq,conseq,timseq

    
    
################################################################################################################################################################
#
#   LOAD VISUAL STUFF
#
################################################################################################################################################################


def loadvisuals(win,cwd,masksize,stimsize,crosssize,factor,xcoo,ycoo,xyconfigs,trial,pos=4):
    #Text
    instructions = visual.TextStim(
        win=win,
        units='norm',
        height = 0.1,
        wrapWidth=1.7,
        alignHoriz='center', 
        text=  "Press space to begin")
    feedback = visual.TextStim(win=win, units='norm', alignHoriz='center', height = 0.1, text= "End of block ")
    infoscr = visual.TextStim(win=win,units='norm', height = 0.1, wrapWidth=1.7,
        alignHoriz='center', text=  "[C]alibration\n[B]uttons\n[T]rial number\n[S]tart " + trial[1:].lower())

    # fixation cross
    fixation = visual.ShapeStim(win, 
        vertices=((0, -0.5*crosssize), (0, 0.5*crosssize), (0,0), (-0.5*crosssize,0), (0.5*crosssize, 0)),
        lineWidth=round(crosssize/8),
        closeShape=False,
        lineColor="red")
    
    # fixation cross
    fixation2 = visual.ShapeStim(win, 
        vertices=((0, -0.5*crosssize), (0, 0.5*crosssize), (0,0), (-0.5*crosssize,0), (0.5*crosssize, 0)),
        lineWidth=round(crosssize/8),
        closeShape=False,
        lineColor="white")    
    ###############################################################################    
    # prepare cues
    r=crosssize/2.5
    cuecolors = [[1,-1,-1],[-1, -0.06, 1],[-1,-1,-1]] # test for equal luminance Blue [-1,-1,1] 6
    colorsq = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0],[0,1,1,0],[1,0,0,1],[0,0,0,0]]

    ss=10
    xycoo0 = [[xcoo[0]/ss,ycoo[0]/ss],[xcoo[1]/ss,ycoo[1]/ss],[xcoo[2]/ss,ycoo[1]/ss],[xcoo[3]/ss,ycoo[0]/ss]]
    xycoo1 = [[xcoo[0]/ss,ycoo[3]/ss],[xcoo[1]/ss,ycoo[2]/ss],[xcoo[2]/ss,ycoo[2]/ss],[xcoo[3]/ss,ycoo[3]/ss]]
    
    clb0,clb1 = [[visual.BufferImageStim(win)]]*7,[[visual.BufferImageStim(win)]]*7
    
    for c in range(len(clb0)):
        cuelist0=[visual.Circle(win=win)]*4+[fixation]
        for d in range(4):
            cuelist0[d] = visual.Circle(win=win,units='pix',radius=r,fillColor=cuecolors[colorsq[c][d]],lineColor=cuecolors[colorsq[c][d]],pos=(xycoo0[d]))
        clb0[c] = visual.BufferImageStim(win,stim=cuelist0)

    for c in range(len(clb0)):
        cuelist1=[visual.Circle(win=win)]*4+[fixation]
        for d in range(4):
            cuelist1[d] = visual.Circle(win=win,units='pix',radius=r,fillColor=cuecolors[colorsq[c][d]],lineColor=cuecolors[colorsq[c][d]],pos=(xycoo1[d]))
        clb1[c] = visual.BufferImageStim(win,stim=cuelist1)

    cue_list = [clb0,clb1]
    
    ###############################################################################
    # load masks, stimuli, probe
    mask_imgs = []
    for i in os.listdir(cwd+"/masks"):
        if i.startswith ("Mask") and i.endswith(".bmp"):
            mask_imgs.append(cwd+"/masks/"+i)
    
    # create visuals
    mask_display,stim_display = [],[]
    for i in range(pos):
        m = visual.ImageStim(win=win, units = 'pix', size = masksize*factor[i])
        s = visual.TextStim(win=win, units = 'pix', text = '', color= 'red', height = stimsize*1.316*factor[i],font='Arial',bold=False)
        mask_display.append(m) 
        stim_display.append(s)
    stim_display.append(fixation)
    mask_display.append(fixation)
    
    maskcollection=[[],[]]
    for c in range(2):
        for m in range(12):
            shuffle (mask_imgs)
            [mask_display[i].setImage(mask_imgs[i]) for i in range(len(xcoo))]
            [mask_display[i].setPos((xcoo[i],ycoo[xyconfigs[c][i]])) for i in range(len(xcoo))]   
            maskcollection[c].append(visual.BufferImageStim(win,stim=mask_display))
    
    probe_display = visual.TextStim(win=win, units = 'pix', color='red', height=stimsize*1.316*factor[-1], pos=(0,0),font='Arial',bold=False)

#    return instructions,feedback,infoscr,fixation,fixation2,cue_list,maskcollection,stim_display,probe_display
    return instructions,feedback,infoscr,fixation,fixation2,cue_list,mask_imgs,mask_display,stim_display,probe_display

################################################################################################################################################################
#
#   PRESENT STARTING SCREEN
#
################################################################################################################################################################


def infoscreen(win,instructions,infoscrtxt,E1,E2,xcoo,ycoo,crosssize,responsekeys,eye=False,fsample='',fevent=''):        
    wait,testbuttons,getcurtrial = True,False,False
    tn,tntxt="","Enter last trial number: "
    current_trial=0

    infoscrtxt.draw()
    win.flip()
    
    while wait:
        if testbuttons==False and getcurtrial==False:
            for key in event.getKeys():
                if key in ['escape']:
                    win.close()
                    core.quit()
                elif key in ['c']:
                    # validation part
                    win.color=[-0.5]*3
                    win.flip()
                    instructions.setText("Please look at the white dot")
                    instructions.draw()
                    win.flip()
                    calibpos = [[0,0],[-E2,E1],[E2,E1],[-E2,-E1],[E2,-E1],[-E2,0],[0,E1],[0,-E1],[E2,0],
                                [xcoo[0],ycoo[0]],[xcoo[1],ycoo[1]],[xcoo[2],ycoo[1]],[xcoo[3],ycoo[0]],
                                [xcoo[3],ycoo[3]],[xcoo[2],ycoo[2]],[xcoo[1],ycoo[2]],[xcoo[0],ycoo[3]],[0,0]]
                    dot = visual.Circle(win=win, radius=0.15*crosssize, fillColor='white',pos=calibpos[0])
                    core.wait(3)
                    
                    FR = win.getActualFrameRate()
                    if eye: iViewXAPI.iV_ContinueEyetracking()
                    for i in range( len(calibpos)):
                        dot.setPos(calibpos[i])
                        if eye: fsample.write('StartCalibration {ts} Pos {x}\n'.format(ts=i,x=calibpos[i]))
                        for frameN in range (int(FR*3)):
                            dot.draw()
                            win.flip()
                        if eye: fsample.write('EndCalibration {ts}\n'.format(ts = i))
                    if eye: iViewXAPI.iV_PauseEyetracking()
                    win.color=[-1]*3
                    win.flip()
                    infoscrtxt.draw()
                    win.flip()
            
                elif key in ['b']:
                    # test buttons part
                    testbuttons = True
                    instructions.setText("Please press the buttons for 'Present','Not present','Don't know'")
                    instructions.draw()
                    win.flip()
                elif key in ['t']:
                    getcurtrial = True
                    
                    instructions.setText(tntxt+tn)
                    instructions.draw()
                    win.flip()
                elif key in ['s','space','return']:
                    wait = False
        
        elif testbuttons:
            buttontext = ["Present","Not present","Don't know"]
            for key in event.getKeys():
                if key in responsekeys:
                    instructions.setText(buttontext[responsekeys.index(key)%3])
                    instructions.draw()
                    win.flip()
                elif key in ['return']:
                    testbuttons = False
                    infoscrtxt.draw()
                    win.flip()
                elif key in ['escape']:
                    win.close()
                    core.quit()
                else:
                    instructions.setText('what')
                    instructions.draw()
                    win.flip()   
                    print key
        elif getcurtrial:
            for key in event.getKeys():
                if key in ['return']: 
                    current_trial = int(tn)
                    getcurtrial=False
                    infoscrtxt.draw()
                    win.flip()
                elif key in ['delete','backspace']:
                    tn = tn[:-1]
                    instructions.setText(tntxt+tn)
                    instructions.draw()
                    win.flip()
                elif key in ['escape']:
                    win.close()
                    core.quit()
                else:
                    tn = tn+key
                    instructions.setText(tntxt+tn)
                    instructions.draw()
                    win.flip()
    
    return current_trial


















