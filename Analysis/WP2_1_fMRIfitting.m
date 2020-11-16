%=======================================================================================
% WP2_1_fMRIfitting
%=======================================================================================
% PROJECT:      WP2 - VF and eccentricity modulations 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Prepare all fMRI data across all participants for SPM
% OUTPUT:       fMRI data organised in /Data_proc_fun
% -------------------------------------------------------------------------
% 2020.06.23 New
% 190412 Moved from analyse functional
% includes STC
% with explicit mask
% 190612 current ver
%        only model specification and estimation, 
% ---------- Use:
% maindir = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity';
% partdirs = dir([maindir '/Data/TE*/session2/functiona*']);partdirs = {partdirs.folder}';
% partid = cellfun(@(x) x(58:61),partdirs,'UniformOutput',false);
% outdir = cellfun(@(x) [maindir '/GroupResults/Functional/Normalize_3vs2.5mm/' x],partid,'UniformOutput',false);

% ----------
function dirs = WP2_1_fMRIfitting(dirs)
if nargin==6, contrastbool = true; else, contrastbool = false; end
%%
for p = 1:length(dirs.fun.prepped)
% Set variables
rpdir = dir([dirs.fun.prepped{p} '/r*/rp_a4D.txt']); 
swadir = dir([dirs.fun.prepped{p} '/r*/swa4D.nii']); 
condname = 'cond_prb';
conddir = dir([dirs.fun.prepped{p} '/r*/' condname '.mat']);
outdir = [dirs.fun.prepped{p} '/out'];
if exist(outdir,'dir')==0,mkdir(outdir);else,warning('%s already exists!',outdir);end

runs=1:length(swadir); 
vols=170;
swa_runs = cell(size(runs));
for r = runs
    for v = 1:vols
        swa_runs{r}{v} = [swadir(r).folder filesep swadir(r).name ',' num2str(v)];
    end
end
disp(cellfun(@length,swa_runs))

tname = '';
tconfunc = '';

% ---------------------------------% ---------------------------------
% MODEL SPECIFICATION AND ESTIMATION
spm('defaults','fMRI');
if exist([outdir '/SPM.mat'],'file') == 0
disp('Performing model specification and estimation...');
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 52;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 40; %second slice in find([0:2000/26:1999 0:2000/26:1999]==1000) from STC
for r = runs
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = swa_runs{r}';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {[conddir(r).folder filesep conddir(r).name]};
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {[rpdir(r).folder filesep rpdir(r).name]};
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
end
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {explmask};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
spm_jobman('run',matlabbatch); 

% ---------------------------------% ---------------------------------
% ESTIMATE

clear matlabbatch
spm_figure('GetWin','Graphics');
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[outdir '/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch); 
disp('Model specification and estimation finished!');
disp('Check and save spm.ps file in Matlab folder');
clear r;

% disp('Check SPM.mat')
% load([outdir '/SPM.mat'])
% % SPM.xY.P -> input files
% % SPM.Sess.U -> regressors
% fprintf('Runs: %d (%s)\nTotal volumes: %d\n%s\n\n',length(SPM.nscan),num2str(SPM.nscan),size(SPM.xY.P,1),SPM.xY.P(1,:))
% disp(SPM.xX.name')

% ---------------------------------% ---------------------------------
% CONTRAST 

if contrastbool
    
    load([outdir '/SPM.mat'])
    if isfield(SPM,'xCon')
        warning('SPM.mat already holds estimated contrasts!')
    else
        disp('Performing F contrast estimation...')
        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = {[outdir '/SPM.mat']};
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'F_all';
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = diag(contains(SPM.xX.name,'bf(1)')*1); % +(1) for only hrf
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 0;
        spm_jobman('run',matlabbatch); 

        if ~isempty(tconfunc)
            tcon = cellfun(@(x) x(SPM.xX.name),tconfunc,'UniformOutput',false);

            disp('Performing T contrast estimation...')
            clear matlabbatch
            matlabbatch{1}.spm.stats.con.spmmat = {[outdir '/SPM.mat']};
            for i=1:length(tname)
                matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = tname{i};
                matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = tcon{i};
                matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
            end
            matlabbatch{1}.spm.stats.con.delete = 0;
            spm_jobman('run',matlabbatch); 
            clear i;
        end
        
        disp('Contrast estimation finished!')
    end 

end
end
%% 200819 Add interaction between eccentricity and #VF

load([outdir '/SPM.mat'])
if length(SPM.xCon) == 5
tname = {'L5R5B10>L10R10B5'};
tcon = {contains(SPM.xX.name,strsplit('L5*bf(1) R5*bf(1) B10*bf(1)')) - ...
    contains(SPM.xX.name,strsplit('L10*bf(1) R10*bf(1) B5*bf(1)'))};

disp('Performing T contrast estimation...')
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = {[outdir '/SPM.mat']};
for i=1:length(tname)
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = tname{i};
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = tcon{i};
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
end
matlabbatch{1}.spm.stats.con.delete = 0;
spm_jobman('run',matlabbatch); 
disp('T contrast estimation of #VF x Ecc Interaction finished!')
end
%dirs.fun.con{end+1} = conddir;
end