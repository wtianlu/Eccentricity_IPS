% WP2_2_analyse_ROI
% 2021.05.08 Revision
% Copied from: WP2_2_analysefMRI.m 
% Get ROIs: line 34 onwards
% Check PSTH results: line 249 onwards

%% Get ROIs from Kastner's max probability maps

spm('defaults','FMRI');

tplfn = {'/home/tianlu/Downloads/ProbAtlas_v4/subj_vol_all/maxprob_vol_lh.nii';
    '/home/tianlu/Downloads/ProbAtlas_v4/subj_vol_all/maxprob_vol_rh.nii'};
outdir = [dirs.fun.main '/Brain_VOI_Wang_Kastner/210509_IPSV1'];

if exist(outdir,'dir')==0,mkdir(outdir);end

roilab = {'[18,19]','[21,22]','[1,2]'}; %IPS0/1, IPS3/4, V1v/V1d
roiname = {'IPS01','IPS34','V1'};
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
%% Get ROIs from Kastner's max probability maps method 2 % line 73

fn = '/home/tianlu/Downloads/ProbAtlas_v4/subj_vol_all/perc_VTPM_vol_roi1_lh.nii';

V = spm_vol(fn);
data = spm_read_vols(V);

% Get values
peak90val = find(data>round(max(data(:))*0.9)); % select values over 90% of the max value
allval = data(data>0); allval = sort(allval,'descend');
peak90all = allval(1:round(length(allval)*0.1)); % select values in the top 90% of all values
peakvals = (allval>round(max(allval)*0.9));

hist(allval);hold on
[n,h] = hist(allval(allval>round(max(allval)*0.9)));
bar(h,n,'red')



%% Check PSTH results
clc; close all
% load('/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_fun/ROI/smoothed_noint/psthresult200810.mat')
load('/home/tianlu/Documents/Projects/2018_fMRI_ECC/Data_proc_fun/ROI/smoothed_noint/psthresult200810.mat')

roinames = strsplit('IPS01lh IPS01rh IPS34lh IPS34rh'); 
% 3/4 = middle IPS; 0/1 are the posterior IPS 


psthdata = cell(np,4,6); % p, r, c
for p = 1:np
    for r = 1:4
        for c = 1:6 
            data = psthresult{p,r}(c,:);
            psthdata{p,r,c} = mean(horzcat(data{:}),2);
        end
    end
end


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
% legend(condtitles([2 3 6 1 4 5]));

%% === Calculate anova ===

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
