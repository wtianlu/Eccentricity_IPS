% WP2_2_analyseFunc
%=======================================================================================
% PROJECT: VF and eccentricity modulations 
% AUTHOR: Lulu Wang
% INSTITUTION: KU Leuven
% AIM: Analyse beta maps and SPM.mat
% =======================================================================================
% 2020.08.10 middle and posterior IPS, check unilateral vs bilateral 

% Data plan (2020.05.30)
% 1. Get IPS0/1, FEF, and SPL from Kastner's map, convert to .mat using
%    SPM result
% 2. Extract PSTH from ROIs with non-separated cues (6 conditions)

%%
function dirs = WP2_2_analysefMRI(dirs)

%% Initialisation
fprintf('\n\n***************** Start WP2_2_analysefMRI *****************\n\n')

% ----------------------------------------------------------------------------------------
% Load participant data
pinfo = readtable([dirs.beh.main '/participant_info.xlsx']);
p_incl = pinfo.IDCode(~contains(pinfo.Other,'Excluded:')); 
np = length(p_incl);


%% Get ROIs from Kastner's max probability maps

spm('defaults','FMRI');

tplfn = {'/Users/Lulu/Documents/Documents:Software/Brain VOI Wang_Kastner/maxprob_vol_lh.nii';
'/Users/Lulu/Documents/Documents:Software/Brain VOI Wang_Kastner/maxprob_vol_rh.nii'};
outdir = '/Users/Lulu/Documents/Documents:Software/Brain VOI Wang_Kastner/200706_IPS01-34';
if exist(outdir,'dir')==0,mkdir(outdir);end

roilab = {'[18,19]','[21,22]'}; %IPS0/1, IPS3/4
roiname = {'IPS01','IPS34'};
hemname = cellfun(@(x) {x(end-5:end-4)},tplfn);

for h = 1:length(tplfn)
    for r = 1:length(roilab)
        roiout = sprintf('%s/%s%s_mni',outdir,roiname{r},hemname{h});
        
        if exist([roiout '.nii'],'file')==0
            clear matlabbatch
            matlabbatch{1}.spm.util.imcalc.input = tplfn(h);
            matlabbatch{1}.spm.util.imcalc.output = roiout;
            matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
            matlabbatch{1}.spm.util.imcalc.expression = sprintf('ismember(i1,%s)',roilab{r});
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = -7;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run',matlabbatch); 
        else
            fprintf('%s.nii already created!\n',roiout)
        end
    end
end

%% Get ROIs from full probability maps

fn = '/Users/Lulu/Documents/Documents:Software/Brain VOI Wang_Kastner/perc_VTPM_vol_roi18_lh.nii';
V = spm_vol(fn);
data = spm_read_vols(V);

% Get values
peak90val = find(data>round(max(data(:))*0.9));
allval = data(data>0); allval = sort(allval,'descend');
peak90all = allval(1:round(length(allval)*0.1));
peakvals = (allval>round(max(allval)*0.9));

hist(allval);hold on
[n,h] = hist(allval(allval>round(max(allval)*0.9)));
bar(h,n,'red')

% Using SPM's GUI: 
% /Users/Lulu/Documents/Documents:Software/Brain VOI Wang_Kastner/200706_IPS01-34/perc/imCalc_IPS01lh_mni.mat

%% Check contrasts for second level analysis

for p = 1:length(p_incl)
    clear SPM
    fnspm = sprintf('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/%s/out/SPM.mat',p_incl{p});
    load(fnspm)
    if 0
    if strcmp(SPM.xCon(2).name,'L10R10B10>L5R5B5')
        fprintf('%s contrast 2 eccentriciy: OK\n',p_incl{p})
    else
        fprintf('%s contrast 2 eccentriciy: missing!!!\n',p_incl{p});pause
    end
    else
        if strcmp(SPM.xCon(4).name,'B5B10>L5L10R5R10')
        fprintf('%s contrast 4 uni/bilateral: OK\n',p_incl{p})
    else
        fprintf('%s contrast 4 uni/bilateral: missing!!!\n',p_incl{p});pause
        end
    end
end
%% Calculate second-level 

maindir = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/';
connames = {'Ecc','VFlat','EccXVF'};
confiles = {cellfun(@(x) {[maindir 'Data_proc_fun/' x '/out/con_0002.nii']},p_incl);
            cellfun(@(x) {[maindir 'Data_proc_fun/' x '/out/con_0004.nii']},p_incl);
            cellfun(@(x) {[maindir 'Data_proc_fun/' x '/out/con_0006.nii']},p_incl)};
        
for i = 1:length(connames)
    outdir = ['/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/WholeBrain/' connames{i}];
    if isdir(outdir)==0
        mkdir(outdir);
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = confiles{i};
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outdir '/SPM.mat']};
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat = {[outdir '/SPM.mat']};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'pos';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = '1';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'neg';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = '-1';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess = {};
        matlabbatch{3}.spm.stats.con.delete = 0;
        spm_jobman('run',matlabbatch); 
    end
end


%% Get peak psth current august
clc; 
% fn_roi = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/90%thresholded/IPS0134.mat';
% outfile = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/90%thresholded/psthresult200810.mat';
fn_roi = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/smoothed_noint/IPS0134.mat';
outfile = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/smoothed_noint/psthresult200810.mat';

if exist(outfile,'file')==0
    tic;
    addpath('/Users/Lulu/Documents/MATLAB/Céline SCRIPTS')
    roilabels = {'IPS01','IPS34'};
    outdir = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI';
    if isdir(outdir)==0,mkdir(outdir);end
    load(fn_roi);
    allrois = {IPS01lh,IPS01rh,IPS34lh,IPS34rh};

    % Calculate psth
    psthresult = cell(np,4);
    for p = 1:np
        for r = 1:4
            fprintf('Calc_PSR %s roi %d\n',p_incl{p},r)
            load(sprintf('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/%s/out/SPM.mat',p_incl{p}))
            psthresult{p,r} = CN5_calc_PSR_celine(allrois{r},SPM);
        end
    end
    t2 = toc; fprintf('Elapsed time: %dm %ds\n',floor(t2/60),round(mod(t2,60)))
    save(outfile,'psthresult')
end
rmpath('/Users/Lulu/Documents/MATLAB/Céline SCRIPTS')

%% imclose
if exist('allrois','var')==0
load('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/90%thresholded/IPS0134.mat');
allrois = {IPS01lh,IPS01rh,IPS34lh,IPS34rh};
end

smoothedrois = cell(1,4);
Vref = spm_vol('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/90%thresholded/IPS01lh.nii');
roinames = strsplit('IPS01lh IPS01rh IPS34lh IPS34rh');

close all;
figure('pos',[1         368        1440         437])
sp = [3 5];
for r = 1:4
    rawimage = zeros(63,76,63);
    for i = 1:size(allrois{r},2)
        rawimage(allrois{r}(1,i),allrois{r}(2,i),allrois{r}(3,i)) = 1;
    end
    se = strel('disk',5); %se2 = strel('disk',1);
    climage = imclose(rawimage,se);%imclose(imclose(rawimage,se),se2);
    
    cl2image = imclose(climage,se);
    
    % Check results
    zslices = [35    39    43    47    51];
    for zslice = 1:length(zslices)
        subplot(sp(1),sp(2),zslice)
        imshow(rawimage(:,:,zslices(zslice)))
        subplot(sp(1),sp(2),zslice+5)
        imshow(climage(:,:,zslices(zslice)))
        subplot(sp(1),sp(2),zslice+10)
        imshow(cl2image(:,:,zslices(zslice)))
    end
    
    [i1,i2,i3] = ind2sub(size(rawimage),find(climage));
    smoothedrois{r} = [i1 i2 i3]';
    
    Vnew = Vref;
    Vnew.fname = strrep(strrep(Vref.fname,'90%thresholded','smoothed'),'IPS01lh',roinames{r});
    spm_write_vol(Vnew,climage);
    
end
save('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/smoothed/IPS0134.mat','smoothedrois');

% remove overlapping voxels
lhint = intersect(smoothedrois{1}',smoothedrois{3}','rows');
rhint = intersect(smoothedrois{2}',smoothedrois{4}','rows');
scheckrois = {setdiff(smoothedrois{1}',lhint,'rows')',setdiff(smoothedrois{2}',rhint,'rows')',...
    setdiff(smoothedrois{3}',lhint,'rows')',setdiff(smoothedrois{4}',rhint,'rows')'};
for r = 1:4
    rawimage = zeros(63,76,63);
    for i = 1:size(scheckrois{r},2)
        rawimage(scheckrois{r}(1,i),scheckrois{r}(2,i),scheckrois{r}(3,i)) = 1;
    end
%     se = strel('disk',2);
%     climage = imclose(rawimage,se);
    Vnew = Vref;
    Vnew.fname = strrep(strrep(Vref.fname,'90%thresholded','smoothed_noint'),'IPS01lh',roinames{r});
    spm_write_vol(Vnew,climage);
    
end
save('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/smoothed_noint/IPS0134.mat','scheckrois');


%% Check PSTH results
clc; close all
load('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/smoothed_noint/psthresult200810.mat')

psthdata = cell(np,4,6); % p, r, c
for p = 1:np
    for r = 1:4
        for c = 1:6
            data = psthresult{p,r}(c,:);
            psthdata{p,r,c} = mean(horzcat(data{:}),2);
        end
    end
end

% % Visualise time course - figure per cue
% figure('pos',[357         249        1084         549]) % for slimmer:  327         391        1114         407
% condtitles = strsplit('L10 L5 R5 R10 B10 B5');
% roinames = {'left pIPS','right pIPS','left mIPS','right mIPS'};
% roiplotlab = {'^--','Color',[0 114 178]/255,'LineWidth',1.5; % pIPSlh pIPSrh mIPSlh mIPSrh
%     'o:','Color',[86 180 233]/255,'LineWidth',1.5;
%     '^--','Color',[213 94 0]/255,'LineWidth',1.5;
%     'o:','Color',[230 159 0]/255,'LineWidth',1.5}; 
% spseq = [4 1 2 5 6 3];
% for c = 1:6
%     subplot(2,3,spseq(c));hold on
%     for r = 1:4
%         data = horzcat(psthdata{:,r,c});
%         e = errorbar(1:16,mean(data,2),std(data,[],2),roiplotlab{r,:});
%         e.MarkerFaceColor = e.Color;
%     end
%     xlim([.5 16.5])
%     ylim([-0.15 0.45])
%     xlabel('time from onset (s)'); ylabel('BOLD % signal change')
%     title(condtitles{c})
%     if c==5, lg = legend(roinames); lg.Orientation = 'horizontal';end
% end

% Figure per region
figure('pos',[440   368   723   430])
cueplotlab = {'^--','Color',[86 180 233]/255,'LineWidth',1.5; % l10 l5 r5 r10 b10 b5
    'o:','Color',[86 180 233]/255,'LineWidth',1.5;
    'o:','Color',[213 94 0]/255,'LineWidth',1.5;
    '^--','Color',[213 94 0]/255,'LineWidth',1.5;
    '^--','Color',[155 197 61]/255,'LineWidth',1.5;
    'o:','Color',[155 197 61]/255,'LineWidth',1.5}; 
for r = 1:4
    subplot(2,2,r); hold on
    for c = [2 3 6 1 4 5]
        data = horzcat(psthdata{:,r,c});
        e = errorbar(1:16,mean(data,2),std(data,[],2),cueplotlab{c,:});
        e.MarkerFaceColor = e.Color;
    end
    xlim([.5 16.5]); ylim([-0.15 0.45])
    xlabel('time from onset (s)'); ylabel('BOLD % signal change')
    title(roinames{r})
end
legend(condtitles([2 3 6 1 4 5]));

% === Calculate anova ===

maxpsth = cellfun(@max,psthdata); % p, roi, c
% per ROI separately per hemifield
factors = table(categorical([1 1 2 2 3 3])',categorical([1 0 0 1 1 0])','VariableNames',strsplit('vf ecc'));
for r = 1:4
    t = array2table(squeeze(maxpsth(:,r,:)));
    fprintf('ROI %d %s\n',r,roinames{r})
    rm = fitrm(t,'Var1-Var6~1','WithinDesign',factors);
    ranovatbl = ranova(rm,'WithinModel','ecc*vf');
    etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
    ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
    disp(ranovatbl(:,[1:5 end]))
%     mc = multcompare(rm,'vf','by','ecc');%,'ComparisonType','bonferroni');
mc = multcompare(rm,'vf');    
disp(mc(mc.Difference>0,:))
end

fprintf('\n\nPer roi with hemifield as factor\n')
factors = table(categorical([1 1 2 2 3 3 1 1 2 2 3 3])',categorical([1 0 0 1 1 0 1 0 0 1 1 0])',...
    categorical([0 0 0 0 0 0 1 1 1 1 1 1])', 'VariableNames',strsplit('vf ecc hem')); %left 0 right 1
for r = [1 3]
    t = array2table([squeeze(maxpsth(:,r,:)) squeeze(maxpsth(:,r+1,:))]);
    fprintf('ROI %d %s\n',r,roinames{r})
    rm = fitrm(t,'Var1-Var12~1','WithinDesign',factors);
    ranovatbl = ranova(rm,'WithinModel','ecc*vf*hem');
    etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
    ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
    disp(ranovatbl(:,[1:5 end]))
%     mc = multcompare(rm,'vf','by','ecc');%,'ComparisonType','bonferroni');
%     disp(mc(mc.Difference>0,:))
end

% % per eccentricity
% factors = table(categorical([1 2 3 1 2 3])', categorical([0 0 0 1 1 1])',...
%     'VariableNames',strsplit('vf hem')); %left 0 right 1
% eccols = {[2 3 6],[1 4 5]};
% for i = 1:2
%      t = array2table([squeeze(maxpsth(:,r,eccols{i})) squeeze(maxpsth(:,r+1,eccols{i}))]);
%     fprintf('Ecc %d\n',i)
%     rm = fitrm(t,'Var1-Var6~1','WithinDesign',factors);
%     ranovatbl = ranova(rm,'WithinModel','vf*hem');
%     etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
%     ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
%     disp(ranovatbl(:,[1:5 end]))   
%         mc = multcompare(rm,'vf');%,'ComparisonType','bonferroni');
%     disp(mc(mc.Difference>0,:))
% end

%% Examine LIdata as function of cue - per ROI (1)
% L10 L5 R5 R10 B10 B5
clc;close all;
for r = 1
fprintf('ROI %d\n',r)    
tmp = cellfun(@(x) {x(1:6)}, LIdata(:,r));tmp = horzcat(tmp{:})';

rmtab = array2table(tmp(:,[2 1 3 4 6 5]));
factors = table(categorical([1 1 2 2 3 3]'),categorical([1 2 1 2 1 2])',...
    'VariableNames',strsplit('VF Ecc'));
rm = fitrm(rmtab,'Var1-Var6~1','WithinDesign',factors);
[ranovatbl] = ranova(rm,'WithinModel','VF*Ecc')

figure('pos',[-440   439   250*1.5   159])
offs = .4; violinxticks = [1 1 3 3 5 5];
plot([0 max(violinxticks)+1],[0 0],':','Color',ones(1,3)*.5);hold on
v = violinplot(tmp,violinxticks+[-offs offs -offs offs -offs offs]);
for vv = 1:2:length(v),v(vv).ViolinColor = ones(1,3)*.2;end
for vv = 2:2:length(v),v(vv).ViolinColor = ones(1,3)*.6;end
ylabel('LI');%ylim([-.5 .5])
xticks(unique(violinxticks));xticklabels(strsplit('LVF RVF BVF'));
box off
end

%% L10 L5 R5 R10 no BVF
clc;close all;
for r = 1
fprintf('ROI %d\n',r)    
tmp = cellfun(@(x) {x(1:4)}, LIdata(:,r));tmp = horzcat(tmp{:})';

rmtab = array2table(tmp(:,[2 1 3 4]));
factors = table(categorical([1 1 2 2]'),categorical([1 2 1 2])',...
    'VariableNames',strsplit('VF Ecc'));
rm = fitrm(rmtab,'Var1-Var4~1','WithinDesign',factors);
[ranovatbl] = ranova(rm,'WithinModel','VF*Ecc')

figure('pos',[-440   439   250   159])
offs = .4; violinxticks = [1 1 3 3];
plot([0 max(violinxticks)+1],[0 0],':','Color',ones(1,3)*.5);hold on
v = violinplot(tmp,violinxticks+[-offs offs -offs offs]);
for vv = 1:2:length(v),v(vv).ViolinColor = ones(1,3)*.2;end
for vv = 2:2:length(v),v(vv).ViolinColor = ones(1,3)*.6;end
ylabel('LI');%ylim([-.5 .5])
xticks(unique(violinxticks));xticklabels(strsplit('LVF RVF'));
box off
end





















