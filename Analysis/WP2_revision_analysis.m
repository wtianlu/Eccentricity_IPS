%=======================================================================================
% WP2_revision_analysis
%=======================================================================================
% PROJECT: VF and eccentricity modulations 
% AUTHOR: Lulu Wang
% INSTITUTION: KU Leuven
% CONTENT: Main organisation script for preprocessing, analysis, visualisation
% =======================================================================================
% EDITS:
% 
% 2021.05.08 Recover original results on Ubuntu
%% Initialisation
dirs.main = '/home/tianlu/Documents/Projects/2018_fMRI_ECC';
dirs.beh.main = [dirs.main '/Data_proc_beh'];
dirs.fun.main = [dirs.main '/Data_proc_fun'];
cd([dirs.main '/Code/Eccentricity_IPS/Analysis'])

%% Check behavioural data
% No need to prep behavioural data first, use EccVF_rawdata.txt
dirs = WP2_1_analyse_beh(dirs,qc); 

% Compare results to JASP: perf_s1_cued_no_em.txt, perf_s2_cued_no_em.txt, perf_s1_base_no_em.txt
% Participants without eye tracking data:
% S1: TE19, not enough per condition: TE10 (baseline: TE19/ not enough: TE10,TE25)
% S2: TE08, TE09; not enough per condition: TE20

%% Check neural data
% First get Kastner's ROI for V1 too in order to extract time series
% -> WP2_analyse_ROI.m section 1

% Recreate condition mat files to hold trials with fixation control


dirs = WP2_2_analysefMRI(dirs);




