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

VERSION 190315 - automatically updating filename in alphabet
VERSION 190501 - run script separately per run and new logging
"""
import os
import math
from psychopy import visual, event, core

def getdisplaysettings():
    
    factor = [1.617, 1, 1,1.617,0.3923] # including S and N: 1.65
    sizes = [2.24, 2.77, 0.4, 5, 10]
    dimpx = [round(math.tan(math.radians(d))*620/0.15) for d in sizes] #.18 large (still full) stim, .05 large cue, .15 regular
    stimsize,masksize,crosssize,E1,E2 = dimpx
    responsekeys = ['r','b','y']   
        
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

def loadFnames(subj_id,run,dataFolder):
    if int(subj_id)<10 and len(subj_id)<2:
        subj_id = '0'+subj_id
    fname = "P" + subj_id + "_run" + str(run) + ".log"
    
    if os.path.isfile(dataFolder + fname):
        print "File "+fname+" already exists"
        tmp = 'abcdefghijklmnopqrstuvwxyz'
        num=len([f for f in os.listdir(dataFolder) if f.startswith(fname[0:8])])
        fname = fname[0:8]+tmp[num]+".log"
    print "Created new file "+fname
    return subj_id, fname

def getEDFname(dataFolder,edfFileName):
    if os.path.isfile(dataFolder + edfFileName):
        print "File "+edfFileName+" already exists" 
        tmp = 'abcdefghijklmnopqrstuvwxyz'
        num=len([f for f in os.listdir(dataFolder) if edfFileName[0:7] in f])
        edfFileName = edfFileName[0:7]+tmp[num]+".EDF"
        print "Created new file "+edfFileName    
    return edfFileName


def loadPseqs():
    cueseq = [4,3,0,3,3,4,3,0,4,0,0,4,5,2,2,5,1,1,5,1,2,2,1,5,3,3,4,4,0,0,4,0,4,3,3,0,5,2,5,2,1,1,1,2,2,5,5,1,4,4,3,3,0,3,0,0,0,3,4,4,5,1,2,5,2,2,2,1,5,5,1,1,4,0,4,4,3,3,4,3,0,0,0,3,5,5,5,1,1,2,2,1,2,1,2,5]        
    prbseq = [3,2,4,4,6,4,3,0,6,0,6,1,2,4,2,6,1,6,1,3,6,1,4,4,7,1,7,0,5,7,5,0,3,3,5,3,7,7,5,2,1,5,7,5,2,1,3,0,3,0,4,3,0,6,4,2,6,0,4,6,6,6,6,4,2,4,3,4,2,0,1,1,0,1,2,7,3,3,5,7,0,7,5,5,2,5,1,1,5,0,7,7,5,2,2,7]
    conseq = [1]*24+[0]*24+[1]*24+[0]*24
    return cueseq,prbseq,conseq


def loadEseqs(cwd,subj_id):
    cueseq,prbseq,conseq=[],[],[]
    f = open(cwd + '/sequences/p' + subj_id + 'sequence_f.txt', 'r')
    tmp = f.readlines()
    f.close()  

    for i in range(len(tmp)):
        cueseq.append(int(tmp[i].split('\t')[0]))
        prbseq.append(int(tmp[i].split('\t')[1]))
        conseq.append(int(tmp[i].split('\t')[2]))
    if min(cueseq)==1:
        cueseq = [i-1 for i in cueseq]
    if min(prbseq)==1:
        prbseq = [i-1 for i in prbseq]
    if min(conseq)==1:
        conseq = [i-1 for i in conseq]    

    return cueseq,prbseq,conseq
    
    
    
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
        text=  "Press space to begin",autoLog=True,name='TextStim_instructions')
    feedback = visual.TextStim(win=win, units='norm', alignHoriz='center', height = 0.1, text= "End of block ",autoLog=True,name='TextStim_feedback')
    infoscr = visual.TextStim(win=win,units='norm', height = 0.1, wrapWidth=1.7,
        alignHoriz='center', text=  "[C]alibration\n[B]uttons\n[T]rial number\n[S]tart " + trial[1:].lower(),autoLog=True,name='TextStim_infoscr')

    # fixation cross
    fixation = visual.ShapeStim(win, 
        vertices=((0, -0.5*crosssize), (0, 0.5*crosssize), (0,0), (-0.5*crosssize,0), (0.5*crosssize, 0)),
        lineWidth=round(crosssize/8),
        closeShape=False,
        lineColor="red",
        name='ShapeStim_redcross')
    
    # fixation cross
    fixation2 = visual.ShapeStim(win, 
        vertices=((0, -0.5*crosssize), (0, 0.5*crosssize), (0,0), (-0.5*crosssize,0), (0.5*crosssize, 0)),
        lineWidth=round(crosssize/8),
        closeShape=False,
        lineColor="white",
        name='ShapeStim_whitecross')    
    ###############################################################################    
    # prepare cues
    r=crosssize/2.5
    cuecolors = [[1,-1,-1],[-1, -0.06, 1],[-1,-1,-1]] # test for equal luminance Blue [-1,-1,1] 6
    colorsq = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0],[0,1,1,0],[1,0,0,1]]
  
    ss=10
    xycoo0 = [[xcoo[0]/ss,ycoo[0]/ss],[xcoo[1]/ss,ycoo[1]/ss],[xcoo[2]/ss,ycoo[1]/ss],[xcoo[3]/ss,ycoo[0]/ss]]
    xycoo1 = [[xcoo[0]/ss,ycoo[3]/ss],[xcoo[1]/ss,ycoo[2]/ss],[xcoo[2]/ss,ycoo[2]/ss],[xcoo[3]/ss,ycoo[3]/ss]]
    clb0,clb1 = [[visual.BufferImageStim(win)]]*6,[[visual.BufferImageStim(win)]]*6

    for c in range(len(clb0)):
        cuelist0=[visual.Circle(win=win)]*4+[fixation]
        for d in range(4):
            cuelist0[d] = visual.Circle(win=win,units='pix',radius=r,fillColor=cuecolors[colorsq[c][d]],lineColor=cuecolors[colorsq[c][d]],pos=(xycoo0[d]))
        clb0[c] = visual.BufferImageStim(win,stim=cuelist0,name='BufferImageStim_cue'+str(c)+'conf0',autoLog=False)
    for c in range(len(clb0)):
        cuelist1=[visual.Circle(win=win)]*4+[fixation] #,units='pix',radius=r,autoLog=False
        for d in range(4):
            cuelist1[d] = visual.Circle(win=win,units='pix',radius=r,fillColor=cuecolors[colorsq[c][d]],lineColor=cuecolors[colorsq[c][d]],pos=(xycoo1[d]))
        clb1[c] = visual.BufferImageStim(win,stim=cuelist1,name='BufferImageStim_cue'+str(c)+'conf1',autoLog=False)
        
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
        m = visual.ImageStim(win=win, units = 'pix', size = masksize*factor[i],name='ImageStim_mask'+str(i),autoLog=True)
        s = visual.TextStim(win=win, units = 'pix', text = '', color= 'red', height = stimsize*1.316*factor[i],font='Arial',bold=False,name='TextStim_stim'+str(i),autoLog=True)
        mask_display.append(m) 
        stim_display.append(s)
    stim_display.append(fixation)
    mask_display.append(fixation)
    probe_display = visual.TextStim(win=win, units = 'pix', color='red', height=stimsize*1.316*factor[-1], pos=(0,0),font='Arial',bold=False,name='TextStim_probe',autoLog=True)

    return instructions,feedback,infoscr,fixation,fixation2,cue_list,mask_imgs,mask_display,stim_display,probe_display


################################################################################################################################################################
#
#   PRESENT STARTING SCREEN
#
################################################################################################################################################################


def infoscreen(win,instructions,infoscrtxt,E1,E2,xcoo,ycoo,crosssize,responsekeys,eye=False,current_trial=0):

    if eye: from pylink import *
        
    wait,testbuttons,getcurtrial = True,False,False
    tn,tntxt="","Enter last trial number: "
    
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
                    instructions.setText("Please look at the dot")
                    instructions.draw()
                    win.flip()
                    calibpos = [[0,0],[-E2,E1],[E2,E1],[-E2,-E1],[E2,-E1],[-E2,0],[0,E1],[0,-E1],[E2,0],
                                [xcoo[0],ycoo[0]],[xcoo[1],ycoo[1]],[xcoo[2],ycoo[1]],[xcoo[3],ycoo[0]],
                                [xcoo[3],ycoo[3]],[xcoo[2],ycoo[2]],[xcoo[1],ycoo[2]],[xcoo[0],ycoo[3]],[0,0]]
                    dot = visual.Circle(win=win, radius=0.15*crosssize, fillColor='white',pos=calibpos[0])
                    if eye:
                        getEYELINK().startRecording(1,1,1,1)
                        getEYELINK().sendMessage("CA start")
                    core.wait(3)

                    for i in range( len(calibpos)):
                        if eye: getEYELINK().sendMessage("StartCA "+str(i))
                        dot.setPos(calibpos[i])
                        dot.draw()
                        win.flip()    
                        core.wait(3)
                        if eye: getEYELINK().sendMessage("EndCA "+str(i))
                    if eye:
                        getEYELINK().sendMessage("CA stop")
                        getEYELINK().stopRecording()   
                        getEYELINK().closeDataFile()
                    
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



import numpy as np
import matplotlib.pyplot as plt

def readlogging_run(fname,showplots=False):
#fname = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Paradigm/Session 2/edfpyData/P01_run4.log'
    with open(fname,'r') as f: 
        logdata = f.readlines()   
        
    run = int(fname[-5])
    if run == 0: txtfname = fname[:-8]+'fPractice.txt'
    else: txtfname = fname[:-8]+'fTask.txt'
    if run>1:fwf = 'a+'
    else: fwf = 'w'
    
    #%% Analyse data timings
        
    # Run information in direct logging
    all_ll = [s for s in logdata if 'll_' in s]
    trial_ll = [s for s in all_ll if 'TRIAL' in s]
    print "\n".join(s for s in all_ll if s not in trial_ll)
    
    ll_trialstart = np.array([float(l.split()[0]) for l in trial_ll[0::2]])
    ll_trialend = np.array([float(l.split()[0]) for l in trial_ll[1::2]])
    if showplots:
        plt.hist(ll_trialend-ll_trialstart)
        plt.title('ll_trial duration')
        plt.show()
    
    file=open(txtfname,fwf)
    file.write( "\n".join(s for s in all_ll[0:2]))
    file.close()  
    
    ###############################################################################
    # Run start to end in win flips
    all_wf = [s for s in logdata if 'wf_' in s]
    
    run_wf = [[s for s in all_wf if 'wf_RUN' in s],[s for s in all_wf if 'wf_show_feedback' in s]]
    trial_wf = [[s for s in all_wf if 'show_ShapeStim_redcross 0' in s],[s for s in all_wf if 'show_BufferImageStim_cue' in s],
                [s for s in all_wf if 'show_ShapeStim_redcross 1' in s],[s for s in all_wf if 'show_STIM' in s],
                [s for s in all_wf if 'show_MASK' in s],[s for s in all_wf if 'show_ShapeStim_redcross 2' in s],
                [s for s in all_wf if 'show_TextStim_probe' in s],[s for s in all_wf if 'show_ShapeStim_whitecross - REPORTING' in s]]
    break_wf = [[s for s in all_wf if 'show_white cross start - 12t' in s],[s for s in all_wf if 'show_white cross end - 12t' in s]]
    
    print 'Recorded trial win flips: '+' '.join(str(len(l)) for l in trial_wf)
    print 'Recorded break win flips: '+' '.join(str(len(l)) for l in break_wf)
    
    # Check trial components
    trialtimeswf = np.array([[float(s.split()[0]) for s in col] for col in trial_wf])
    if showplots:
        for i in range(len(trial_wf)-1):
            dur = np.round((trialtimeswf[i+1,:]-trialtimeswf[i,:])*1000,1)
            plt.hist(dur)
            plt.title('{} - M:{:.1f}, sd:{:.1f}, range: {} - {}'.format(i+1,np.mean(dur),np.std(dur),min(dur),max(dur)))
            plt.show()
        # ll_trialend-trialtimeswf[-1,:] # difference is 1 frame
        # ll_trialstart - trialtimeswf[0,:] # difference is also about 1 frame
        dur = (trialtimeswf[0,1:]-trialtimeswf[-1,:-1])*1000 
        plt.plot(dur,marker='o')
        plt.title('{} - M:{:.1f}, sd:{:.1f}, range: {} - {}'.format(i+2,np.mean(dur),np.std(dur),min(dur),max(dur)))
        plt.show()
        dur2=dur[abs(dur-np.mean(dur))<np.std(dur)]
        plt.plot(dur2,marker='o')
        plt.title('{} - M:{:.1f}, sd:{:.1f}, range: {} - {}'.format(i+9,np.mean(dur2),np.std(dur2),min(dur2),max(dur2)))
        plt.show()        
        plt.plot(np.array([float(l.split()[0]) for l in break_wf[1]])-np.array([float(l.split()[0]) for l in break_wf[0]]),marker='o')
        plt.title('Break time')
        plt.show()    
    #%% Get presented stimuli, cues, responses
    stimtrial = [int(l.split()[-1]) for l in logdata if 'TRIAL_start' in l]
    stimcue = [int(l.split()[-1][-6]) for l in trial_wf[1]]
    stimconf = [int(l.split()[-1][-1]) for l in trial_wf[1]]
    stimletters = [[l.split()[-1][1] for l in logdata if 'TextStim_stim'+str(ll)+': text' in l] for ll in range(4)]
    stimletters = map(list,map(None,*stimletters))
    stimprobe =  [l.split()[-1][1] for l in logdata if 'TextStim_probe: text' in l]
    # get location of probe in stimulus display or return -1 if not present
    probeloc = [stimletters[t].index(stimprobe[t]) if stimprobe[t] in stimletters[t] else -1 for t in range(len(stimprobe))]
    validity = []
    for t in range(len(stimprobe)):
        if probeloc[t]<0:v=2
        elif stimcue[t]<4 and stimcue[t]==probeloc[t]: v=1 
        elif stimcue[t]==4 and probeloc[t] in [0,3]:v=1
        elif stimcue[t]==5 and probeloc[t] in [1,2]:v=1
        else: v=0
        validity.append(v)
    
    # get responses
    tstartind = [logdata.index(l) for l in logdata if 'TRIAL_start' in l]
    tendind = [logdata.index(l) for l in logdata if 'TRIAL_end' in l]
    
    responselist,allrtlist = [],[]
    for i in range(len(tstartind)):
        r = [l for l in logdata[tstartind[i]:tendind[i]] if 'Keypress: ' in l]
        responselist.append([k.split()[-1] for k in r])
        allrtlist.append([float(k.split()[0]) for k in r])
    if showplots:    
        plt.hist([len(s) for s in allrtlist])
        plt.title('Number of responses per trial')
        plt.show()
    
    ########## Check responses
    '''
    -1 if no or invalid response given
    0 if incorrect
    1 if correct
    2 if DK
    3 if multiple responses given
    '''
    trialresponse = []
    for t in range(len(responselist)):
        if len(responselist[t])<1: r = -1
        elif len(responselist[t])>1: r = 3
        elif responselist[t][0]=='y': r = 2
        elif responselist[t][0]=='r':
            if probeloc[t]<0: r=0
            else: r=1
        elif responselist[t][0]=='b':
            if probeloc[t]<0: r=1
            else: r=0
        else: 
            print [t]+responselist[t]
            r=-1
        trialresponse.append(r)
    if showplots:
        plt.hist(trialresponse)
        plt.title('Analysed responses')
        plt.show()
    print 'Percentage correct: {:.1f}%'.format(100.0*trialresponse.count(1)/(trialresponse.count(1)+trialresponse.count(0)))
    
    # Get reaction time of correct trials from probe onset
    rtlist = np.array([(allrtlist[t][0]-trialtimeswf[6,t])*1000 if trialresponse[t]==1 else np.nan for t in range(len(allrtlist))])
    if showplots:
        plt.hist(rtlist[np.isnan(rtlist)==False])
        plt.title('RT correct trials')
        plt.show()
    
    
    #%% Write to txt for further analysis
    file=open(txtfname,fwf)
    file.writelines(['{}\t{}\t{}\t12\t{}\tT:{}\tP:{}\t{}\t{}\tR-{}\t[{}]\tTT:[{}]\n'.format(stimtrial[t],stimcue[t],probeloc[t],stimconf[t],\
                     ''.join(stimletters[t]), stimprobe[t],trialtimeswf[0,t],np.round(1000*(trialtimeswf[4,t]-trialtimeswf[3,t]),1),\
                     ''.join(responselist[t]),','.join([str(rr) for rr in allrtlist[t]]),','.join([str(rr) for rr in trialtimeswf[:,t]]))\
                     for t in range(len(stimcue))])
    
    file.write( "\n".join(s for s in all_ll[-2:]))
    file.close()  
    print 'Log saved to txt'














