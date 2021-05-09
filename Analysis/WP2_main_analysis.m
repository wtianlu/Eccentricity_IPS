%=======================================================================================
% WP2_main_analysis
%=======================================================================================
% PROJECT: VF and eccentricity modulations 
% AUTHOR: Lulu Wang
% INSTITUTION: KU Leuven
% CONTENT: Main organisation script for preprocessing, analysis, visualisation
% =======================================================================================
% EDITS:
% 
% 2021.05.08 Recover original results on Ubuntu, re-do with smaller n
% ---
% 2020.08.10 thesis version 
% 2020.06.15 ms version 20200609
% 2020.05.30 Gathering separate scripts
% ms version 200530_eccentricity_ms

%% ----------------------------------------------------------------------------------------
% 0. Initialise + load raw data  
% ----------------------------------------------------------------------------------------
clc; clearvars; close all;
% dirsmat = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Code/dirs.mat';
dirsmat = '/home/tianlu/Documents/Projects/2018_fMRI_ECC/Code/Eccentricity_IPS/Analysis';
if exist(dirsmat,'file')==0
    dirs.main = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity'; 
    dirs.raw.main = [dirs.main '/Data_raw'];
    dirs.beh.main = [dirs.main '/Data_proc_beh'];
    dirs.fun.main = [dirs.main '/Data_proc_fun'];
else
    load(dirsmat)
end
cd([dirs.main '/Code']) 

disp('Ready!')

%% If finished, save dirsmat
save(dirsmat,'dirs')

%% 1. Check behavioural accuracy data --------------------------------------------------------------
% qc: set to 1 for quality check
qc = 1; 

dirs = WP2_1_prepdata_beh(dirs);
dirs = WP2_1_analyse_beh(dirs,qc); 

% 2. Fit TVA data (with adjustments for C and t0
% dirs.tva = [dirs.main '/Data_tva'];
% dirs = WP2_1_tvafitting(dirs);
% dirs = WP2_2_analyse_tvafit(dirs);

save(dirsmat,'dirs')

%% 2. Check neural and TVA data --------------------------------------------------------------

dirs = WP2_1_prepdata_fMRI(dirs);
dirs = WP2_1_fMRIfitting(dirs);
dirs = WP2_2_analysefMRI(dirs);

save(dirsmat,'dirs')




