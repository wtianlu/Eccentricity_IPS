%=======================================================================================
% WP2_1_analyse_beh
%=======================================================================================
% PROJECT:      WP2 - VF and eccentricity modulations 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Analyse all behavioiral data across participants 
% OUTPUT:       
% -------------------------------------------------------------------------
% 2020.08.10 separate cued and distributed (even leave out distributed
% 2020.06.09 added DK check
% 2020.05.30 
% -------------------------------------------------------------------------
function dirs = WP2_1_analyse_beh(dirs,qc)

if nargin<2, qc=false; end
%% Initialisation
fprintf('\n\n***************** Start WP2_1_prepdata_beh *****************\n\n')

% Read participant information
pinfo = readtable([dirs.beh.main '/participant_info.xlsx']);
idx_incl = ~contains(pinfo.Other,'Excluded:');
p_incl = pinfo.IDCode(idx_incl); 
np = length(p_incl);
fprintf('Included participants: %d\n',np);fprintf('%s\t',p_incl{:});
fprintf('\nAge: %.1f +- %.1f, range %d - %d\n', mean(pinfo.Age(idx_incl)), std(pinfo.Age(idx_incl)), ...
    min(pinfo.Age(idx_incl)), max(pinfo.Age(idx_incl)))

% Load data
fprintf('\nLoading raw data...\n')

bd.raw = readtable([dirs.beh.main '/EccVF_rawdata.txt']);
%  -- Process raw data
% Filter out excluded participants
bd.fil = bd.raw(contains(bd.raw.pidcode,p_incl),:);
% Filter out eye movements
% bd.fil_em = bd.fil(bd.fil.eyedata~=1,:); % 1: eye movement detected; NaN: no available eye data
bd.fil_em = bd.fil(bd.fil.eyedata == 0,:); % 20210508 edit -> also remove trials where fixation cannot be confirmed

fprintf('\nReady!\n')
%% Quality check: DK rates of excluded participants
fprintf('\nChecking DK rates\n')
DKr = splitapply(@(x) sum(x==3)/sum(~isnan(x)),...
    bd.raw.Response(bd.raw.session==1),bd.raw.pidnum(bd.raw.session==1));
% disp(round(DKr*100))
fprintf('Outliers: %.2f\t%.2f\n',DKr([1 8]))
fprintf('DK-rate M+-SD: %.2f +- %.2f, range %.2f - %.2f\n',mean(DKr),std(DKr),min(DKr),...
    max(DKr(~ismember(1:length(DKr),[1 8]))))

%% Quality check: eye tracking
clc
data = bd.fil;%(bd.fil.session==s,:);
naneye = splitapply(@(x) mean(isnan(x)),data.eyedata,findgroups(data.pidnum));
fprintf('Missing eye data: %.1f +- %.1f\n',mean(naneye*100),std(naneye*100))

for s = 1:2
fprintf('\nChecking eye tracking data session %d\n',s)
data = bd.fil(bd.fil.session==s,:);
em = splitapply(@(x) sum(x==1)/sum(~isnan(x)),data.eyedata,findgroups(data.pidnum)); em = em(~isnan(em));
fprintf('Eye movements: %.2f +- %.2f\n',mean(em*100),std(em*100))

data = data(data.CueI<6,:);
emc = splitapply(@(x) sum(x==1)/sum(~isnan(x)),data.eyedata,findgroups(data.pidnum,data.CueI));
emc = reshape(emc,6,18)'; emc(isnan(emc)) = 0;
rmtab = array2table(emc(:,1:6));
factors = table(categorical([2 1 1 2 2 1]'),categorical([1 1 2 2 3 3]'),'VariableNames',...
    strsplit('ecc vf'));
rm = fitrm(rmtab,'Var1-Var6~1','WithinDesign',factors);
ranovatbl = ranova(rm,'WithinModel','ecc*vf');
etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
disp(ranovatbl(:,[1:5 end]))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse cued data separate from baseline data (200810)
close all; clc
no_em = true; % if true: use only the trials with valid eye tracking data

% -- Select SESSION --
figure('pos',[ 724   675   182*3   123*4])

for s = 1:2
    fprintf('\nChecking performance data session %d\n',s)
    p_incl = pinfo.IDCode(idx_incl); 
    
    subplot(1,3,s+(s-1))
    if no_em
        data = bd.fil_em(and(bd.fil_em.session==s,bd.fil_em.CueI<6),:); 
    else
        data = bd.fil(and(bd.fil.session==s,bd.fil.CueI<6),:); 
    end
    np = length(unique(bd.fil_em.pidnum(bd.fil_em.session==s)));
    tmp = setdiff(p_incl,unique(data.pidcode));
    fprintf('Removed participant %s\n',tmp{:})
    p_incl = unique(data.pidcode);


    % Get percentage of correct responses per participant and per cue
    perf = zeros(np,6); dk = zeros(np,6);
    pincln = unique(data.pidnum)';
    for p = 1:length(pincln)
        for c = unique(data.CueI)'
            perf(p,c+1) = mean(data.Correct(and(data.pidnum==pincln(p),data.CueI == c))==1); % including DK
    %         perf(p,c+1) = sum(data.Correct(and(data.pidnum==pincln(p),data.CueI == c))==1)./...
    %             sum(data.Correct(and(data.pidnum==pincln(p),data.CueI == c))> -1); % only trials without DK
    %         dk(p,c+1) = mean(data.Correct(and(data.pidnum==pincln(p),data.CueI == c))==-1);
        end
    end

    % 20210508 Eye tracking data fixes - remove participants where not enough
    % trials per condition remain after removing non-fixation trials
    fprintf('Remove participants with insufficient trials:\n%s\n',p_incl{any(isnan(perf),2)})
    % To use for the functional analyses: 
    p_incl = p_incl(~any(isnan(perf),2));
    if s == 2, save('p_incl_s2_no_em.mat','p_incl');end
    
    perf = perf(~any(isnan(perf),2),:);
    fprintf('Remaining participants: %d\n',size(perf,1))



% -- Calculate 2-way ANOVA

factors = table(categorical([1 1 2 2 3 3])',categorical([1 0 0 1 1 0])','VariableNames',strsplit('vf ecc'));
t = array2table(perf*100); %disp([factors table(round(mean(t{:,:})',1))])
rm = fitrm(t,'Var1-Var6~1','WithinDesign',factors);
ranovatbl = ranova(rm,'WithinModel','ecc*vf');
etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
disp(ranovatbl(:,[1:5 end]))
mc = multcompare(rm,'vf','by','ecc');%,'ComparisonType','bonferroni');
disp(mc(mc.Difference>0,:))
mc = multcompare(rm,'ecc','by','vf');%,'ComparisonType','bonferroni');
disp(mc(mc.Difference>0,:))

tmp = [ranovatbl.Properties.RowNames(3:2:7) num2cell([ranovatbl.DF(3:2:7) ...
    ranovatbl.DF(4:2:8)]) table2cell(ranovatbl(3:2:7,[4 5 9]))]';
fprintf('%s\tF_{\\,%d, %d} = %.2f, p = %.3f, \\eta^2_p = %.3f\n',tmp{:,:})


% Visualise in violinplots

offs = .4; violinxticks = [1 1 3 3 5 5];
if 0
    % violinplot
    v = violinplot(t{:,[2 1 3 4 6 5]},violinxticks+[-offs offs -offs offs -offs offs]); %5 and 10
    for vv = [1 3 5],v(vv).ViolinColor = ones(1,3)*.2;end
    for vv = [2 4 6],v(vv).ViolinColor = ones(1,3)*.6;end
    lg = legend([v(1).ScatterPlot,v(2).ScatterPlot],['5' char(176)],['10' char(176)],'location','northeast'); 
else

%         x=t{:,1};
%         SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([0.025  0.975],height(t)-1);      % T-Score
%         CI = mean(x) + ts*SEM;   

    x = ts(2)*std(t{:,[2 3 6]})./sqrt(height(t)); %std(t{:,[2 3 6]})
    errorbar(unique(violinxticks)-ones(1,3)*offs,mean(t{:,[2 3 6]}),x,...
        'ko-','MarkerFaceColor',ones(1,3)*.9); hold on
    x = ts(2)*std(t{:,[1 4 5]})./sqrt(height(t)); %std(t{:,[1 4 5]})
    errorbar(unique(violinxticks)+ones(1,3)*offs,mean(t{:,[1 4 5]}),x,...
        'ko-','MarkerFaceColor',ones(1,3)*.5); 
end
xticks(unique(violinxticks));xticklabels(strsplit('LVF RVF BVF')); ylabel('% correct')
ylim([30 100])
box off; 

% Save table for JASP analysis
t.Properties.VariableNames = strsplit('L10 L5 R5 R10 B10 B5');
writetable(t,sprintf('perf_s%d_cued_no_em.txt',s));

end


% Check baseline trials 200819
%close all; clc



% -- Select SESSION -- baseline trials

fprintf('\nBaseline trials\n')
if no_em
    data = bd.fil_em(bd.fil_em.runnr==13,:); 
else
    data = bd.fil(bd.fil.runnr==13,:);
end
p_incl = pinfo.IDCode(idx_incl);
tmp = setdiff(p_incl,unique(data.pidcode));
fprintf('Removed participant %s\n',tmp{:})
p_incl = unique(data.pidcode);

% Get percentage of correct responses per participant and per cue
perf = zeros(np,5); dk = zeros(np,5);
pincln = unique(data.pidnum)';
for p = 1:length(pincln)
    for c = unique(data.Probeloc)'
        perf(p,c+1) = mean(data.Correct(and(data.pidnum==pincln(p),data.Probeloc == c))==1); % with DK trials
%         perf(p,c+1) = sum(data.Correct(and(data.pidnum==pincln(p),data.Probeloc == c))==1)./...
%                     sum(data.Correct(and(data.pidnum==pincln(p),data.Probeloc == c))> -1);
        
%         dk(p,c+1) = mean(data.Correct(and(data.pidnum==pincln(p),data.Probeloc == c))==-1);
    end
end
perf = perf(:,2:end); dk = dk(:,2:end); % remove NP trials

% 20210508 Eye tracking data fixes - remove participants where not enough
% trials per condition remain after removing non-fixation trials
fprintf('Remove participants with insufficient trials:\n%s\n',p_incl{any(isnan(perf),2)})
perf = perf(~any(isnan(perf),2),:);
fprintf('Remaining participants: %d\n',size(perf,1))

% Calculate 2-way ANOVA
factors = table(categorical([1 1 2 2])',categorical([1 0 0 1])','VariableNames',strsplit('vf ecc'));
t = array2table(perf*100); disp([factors table(round(mean(t{:,:})',1))])
rm = fitrm(t,'Var1-Var4~1','WithinDesign',factors);
ranovatbl = ranova(rm,'WithinModel','ecc*vf');
etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
disp(ranovatbl(:,[1:5 end]))
% % mc = multcompare(rm,'vf','by','ecc');%,'ComparisonType','bonferroni');
% mc = multcompare(rm,'vf');
% disp(mc(mc.Difference>0,:))
% % mc = multcompare(rm,'ecc','by','vf');%,'ComparisonType','bonferroni');
% mc = multcompare(rm,'ecc');
% disp(mc(mc.Difference>0,:))

tmp = [ranovatbl.Properties.RowNames(3:2:7) num2cell([ranovatbl.DF(3:2:7) ...
    ranovatbl.DF(4:2:8)]) table2cell(ranovatbl(3:2:7,[4 5 9]))]';
fprintf('%s\tF_{\\,%d, %d} = %.2f, p = %.3f, \\eta^2_p = %.3f\n',tmp{:,:})

% -- Visualise

offs = .4; violinxticks = [1 1 3 3];
% figure('pos',[ 724   675   182   123])
subplot(1,3,2)
ts = tinv([0.025  0.975],height(t)-1); 
x = ts(2)*std(t{:,[2 3]})./sqrt(height(t)); %std(t{:,[2 3 6]})
errorbar(unique(violinxticks)-ones(1,2)*offs,mean(t{:,[2 3]}),x,...
    'ko-','MarkerFaceColor',ones(1,3)*.9); hold on

x = ts(2)*std(t{:,[1 4]})./sqrt(height(t)); %std(t{:,[1 4 5]})
errorbar(unique(violinxticks)+ones(1,2)*offs,mean(t{:,[1 4]}),x,...
    'ko-','MarkerFaceColor',ones(1,3)*.5); 
lg = legend(['5' char(176)],['10' char(176)],'location','northeast'); 
xticks(unique(violinxticks));xticklabels(strsplit('LVF RVF BVF')); ylabel('% correct')
ylim([30 100])
box off; legend boxoff
lg.Orientation = 'horizontal';

% Save table for JASP analysis
t.Properties.VariableNames = strsplit('L10 L5 R5 R10');
writetable(t,'perf_s1_base_no_em.txt');


%% Extra: check effect of exposure duration in session 1
clc;

data = bd.fil(and(bd.fil.session==1,bd.fil.CueI<6),:); 
tmp = repmat(0:5,1,6);
tmp2 = reshape(repmat(0:5,6,1),6*6,1);

factors = table(repmat(categorical([1 1 2 2 3 3])',6,1),repmat(categorical([1 0 0 1 1 0])',6,1),...
    categorical(tmp2), 'VariableNames',strsplit('vf ecc exp'));

rm = zeros(18,height(factors));
for p = 1:18
    pdata = data(strcmp(data.pidcode,p_incl{p}),:);
    for c = 1:height(factors)
        tmpdata = pdata(and(pdata.CueI==tmp(c),str2double(pdata.TimI)==tmp2(c)),:);
        rm(p,c) = mean(tmpdata.Correct==1);
    end
end

t = array2table(rm);
rm = fitrm(t,'rm1-rm36~1','WithinDesign',factors);
ranovatbl = ranova(rm,'WithinModel','ecc*vf*exp');
etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
disp(ranovatbl(:,[1:5 end]))

%% check for response differences for trial types v/i/np
clc
for s = 1:2
    data = bd.fil(and(bd.fil.session==s,bd.fil.CueI<6),:); 
    resp = reshape(splitapply(@(x) sum(x==2)/length(x),data.Response,findgroups(data.pidcode,data.Validity)),3,18)';
    [H,P,CI,STATS]=ttest(resp(:,2),resp(:,3))
end

%%
% Save for R analysis:
% ecc = repmat(factors.Ecc',height(t),1); pid = repmat([1:height(t)]',1,8);
% vf = repmat(factors.VF',height(t),1); cond = repmat(factors.Cond',height(t),1);
% perfR = table(pid(:),t{:,:}(:),ecc(:),vf(:),cond(:),'VariableNames',strsplit('pid perf ecc vf cond'));
% writetable(perfR,'WP2_beh_accuracyR.txt')


fprintf('\n\n***************** End WP2_1_prepdata_beh *****************\n\n')
