%=======================================================================================
% WP2_1_prepdata_fMRI
%=======================================================================================
% PROJECT:      WP2 - VF and eccentricity modulations 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Prepare all fMRI data across all participants for SPM
% OUTPUT:       fMRI data organised in /Data_proc_fun
% -------------------------------------------------------------------------
% 2020.06.23 New
%-----------------------------------------------------------------------
% adapted from SPM_preproc_18....
% 190412 adapted and cleaned up from analyse_functional_190309, separated
% modeling and estimation
% 190506 changed from 190412 -> no STC, and remove segmentation and 
% normalisation of the structural image. Turned script into a function
% 190612 turn into 'official' final version: with STC, with all components
%        clean without figures popping up (but saving them), 2.5mm norm
% ---------------------------------% ---------------------------------

function dirs = WP2_1_prepdata_fMRI(dirs)
% Get files
pinfo = readtable([dirs.beh.main '/participant_info.xlsx']);
p_incl = pinfo.IDCode(~contains(pinfo.Other,'Excluded:')); 
dirspart = dir([dirs.fun.main '/TE*']); dirspart = fullfile({dirspart.folder}',{dirspart.name}');
dirspart = dirspart(contains(dirspart,p_incl));

% Set parameters
funcfiles = dir([dirspart{1} '/r*/4D.nii']);
runs=1:length(funcfiles);
vols = 170;
fwhms = 5; % 8

for p = 1:length(dirspart)
    spfname = '%s/r%d/%s4D.%s'; 
    fprintf(spfname,dirspart{p},1,'','nii');fprintf('\n')

spm('defaults', 'FMRI');
spmdir = spm('Dir');
spm_jobman('initcfg');
% ---------------------------------% ---------------------------------
% SLICE TIMING CORRECTION

disp('Performing slice timing correction...');
for r=runs
    if isdir(dirspart{p})==0,mkdir(dirspart{p});end
    rawfunc = cellfun(@(x) {sprintf(spfname,dirspart{p},r,'',['nii,' num2str(x)])},num2cell(1:vols)');
    if exist(sprintf(spfname,dirspart{p},r,'a','nii'),'file') == 0
    clear matlabbatch
    matlabbatch{1}.spm.temporal.st.scans = rawfunc;
    matlabbatch{1}.spm.temporal.st.nslices = 52;
    matlabbatch{1}.spm.temporal.st.tr = 2;
    matlabbatch{1}.spm.temporal.st.ta = 2-2/26; 
    % (TR-TR/slices) Manual p20: If the next two items are entered in ms, this entry will not be used (can be set to 0)
    matlabbatch{1}.spm.temporal.st.so = [0:2000/26:1999 0:2000/26:1999];
    matlabbatch{1}.spm.temporal.st.refslice = 1000; % middle of acquisition
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    spm_jobman('run',matlabbatch);
    else
    fprintf([spfname ' already exists!\n'],dirspart{p},r,'a','nii')
    end
end
fprintf('Finished slice timing correction of the functional images!\n\n');


% ---------------------------------% ---------------------------------
% REALIGNMENT
disp('Performing motion correction of the functional images...'); 
for r=runs
    if exist(sprintf(spfname,dirspart{p},r,'rp_a','txt'),'file')==0
    disp(['Run ' num2str(r) ' motion correction...']); 
    clear matlabbatch  
    afunc = cellfun(@(x) {sprintf(spfname,dirspart{p},r,'a',['nii,' num2str(x)])},num2cell(1:vols)');
    matlabbatch{1}.spm.spatial.realign.estwrite.data = afunc;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2; 
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
    end
    disp(['Finished motion correction of run ' num2str(r)]);
end

% % quality check -> create figures
% clear pi
% sf = 1;
% if exist([partdir '/results'],'dir')==0,mkdir([partdir '/results']);end
% for r=runs
%     fid = fopen([partdir '/functional/run' num2str(r) '/rp_afunctional4D.txt']);
%     rpdata = textscan(fid,[repmat('%f ',[1,6]) ' \n']);
%     fclose(fid);
%     disp(['Run ' num2str(r) ' - Max absolute displacement values: ' num2str(cellfun(@(x) max(abs(x)),rpdata(1:3)),'%0.3f ') ' mm - '...
%     num2str(cellfun(@(x) max(abs(x))*(180/pi),rpdata(4:6)),'%0.3f ') ' degrees'])
%     rplabels = strsplit('x y z pitch roll yaw');
%     figure('pos',[0 0 600 600],'visible','off');
%     subplot(2,1,1)
%     for i=1:3,plot(rpdata{i});hold on;end
%     legend(rplabels(1:3),'Location','southwest')
%     xlabel('Volume');ylabel('mm')
%     title(['Max translation: ' num2str(cellfun(@(x) max(abs(x)),rpdata(1:3)),'%0.3f ')])
%     subplot(2,1,2)
%     for i=4:6,plot(rpdata{i}*(180/pi));hold on;end
%     legend(rplabels(4:6),'Location','southwest')
%     xlabel('Volume');ylabel('degrees')
%     title(['Max rotation: ' num2str(cellfun(@(x) max(abs(x))*(180/pi),rpdata(4:6)),'%0.3f ')])
%     axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 0.98,[partdir(end-12:end-9) ' run ' num2str(r)])
%     if sf        
%         saveas(gcf,[partdir '/results/' partdir(end-12:end-9) '_realign_run-' num2str(r) '.png'])
%     end
% end


% ---------------------------------% ---------------------------------                                     
% SEGMENTATION

fprintf('\n\nPerforming segmentation of the structural image...\n'); 
rawanat = sprintf('%s/session2/structural/structural3D.nii',strrep(dirspart{p},'proc_fun','raw'));
if exist(strrep(rawanat,'structural3D','y_structural3D'),'file')==0
clear matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[rawanat ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmdir '/tpm/TPM.nii,1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmdir '/tpm/TPM.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmdir '/tpm/TPM.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmdir '/tpm/TPM.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmdir '/tpm/TPM.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmdir '/tpm/TPM.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 1; 
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
spm_jobman('run',matlabbatch); 
end
disp('Finished segmentation of the structural image!');


% ---------------------------------% ---------------------------------
% COREGISTRATION
for r=runs
    disp(['Performing spatial coregistration run ' num2str(r) '...']);
    clear matlabbatch
    afunc = cellfun(@(x) {sprintf(spfname,dirspart{p},r,'a',['nii,' num2str(x)])},num2cell(1:vols)');
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {rawanat}; %target, remains stationary
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {sprintf(sprintfname,dirspart{p},r,'meana','nii,1')};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = afunc;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run',matlabbatch);
end
disp('Finished spatial coregistration!');

% ---------------------------------% ---------------------------------
% NORMALIZE

disp('Performing normalisation of the functional images ...');
if exist(sprintf(spfname,dirspart{p},r,'wa','nii'),'file')
    afunc = cellfun(@(y) {cellfun(@(x) {sprintf(sprintfname,dirspart{p},y,'a',...
        ['nii,' num2str(x)])},num2cell(1:vols)')},num2cell(runs));
    clear matlabbatch
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {strrep(rawanat,'structural3D','y_structural3D')};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = vertcat(afunc{:});
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = 2.5*ones(1,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run',matlabbatch); 
end
disp('Finished normalisation of the functional images!');

disp('Performing normalisation of the structural image...')
if exist(strrep(rawanat,'structural3D','wstructural3D'),'file')==0
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {strrep(rawanat,'structural3D','y_structural3D')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[rawanat ',1']};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch); 
end
disp('Finished normalisation of the structural image!');
% ---------------------------------% ---------------------------------
% SMOOTHING

disp('Performing spatial smoothing...');
if exist(sprintf(spfname,dirspart{p},r,'swa','nii'),'file')==0
    afunc = cellfun(@(y) {cellfun(@(x) {sprintf(sprintfname,dirspart{p},y,'wa',...
        ['nii,' num2str(x)])},num2cell(1:vols)')},num2cell(runs));
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = vertcat(afunc{:});
matlabbatch{1}.spm.spatial.smooth.fwhm = ones(1,3)*fwhms;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch); 
end
% ---------------------------------% ---------------------------------
disp('Preprocessing finished!');

dirs.fun.prepped = dirspart;

end





