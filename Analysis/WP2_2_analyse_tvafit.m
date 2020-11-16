% WP2_1_fitTVA
%=======================================================================================
% PROJECT: VF and eccentricity modulations 
% AUTHOR: Lulu Wang
% INSTITUTION: KU Leuven
% AIM: Analyse tva fitted data
% =======================================================================================
% 2020.05.30
% ----------------------------------------------------------------------------------------
function dirs = WP2_2_analyse_tvafit(dirs)

fprintf('\n\n***************** Start WP2_2_analyse_tvafit *****************\n\n')

% Load fit data
fn_fit = [dirs.tva '/fits_summary.xlsx'];
fd_s1 = readtable(fn_fit,'Sheet','S1_fix_t0_C_fin');
fd_s2 = readtable(fn_fit,'Sheet','S2_fix_t0_C_fin');
fd_bl = readtable(fn_fit,'Sheet','s1bl');
allfits = {fd_s1,fd_bl,fd_s2}; % session 1 selective and baseline, session 2 selective

% Load participant data
pinfo = readtable([dirs.beh.main '/participant_info.xlsx']);
p_incl = pinfo.IDCode(~contains(pinfo.Other,'Excluded:')); 
np = length(p_incl);
np = height(fd_s1);

disp('Ready!')

%% Check fit performance
% clc;
checkfits = {fd_s1,fd_bl}; corrRP = zeros(np,2);
for i = 1%:length(checkfits) % Session 2 not discussed in 3.1.3
    fprintf('Session %d\n',i)
    obspredcol = {find(contains(checkfits{1}.Properties.VariableNames,'ObsPr'));
        find(contains(checkfits{1}.Properties.VariableNames,'PredPr'))};
    for p = 1:height(checkfits{i})
        tmp = cellfun(@(x) {checkfits{i}{p,x}},obspredcol);
        tmp = vertcat(tmp{:})';
        [R,P] = corrcoef(tmp);
        corrRP(p,:) = [R(1,2),P(1,2)];
    end
    disp(array2table([mean(corrRP);std(corrRP);min(corrRP);max(corrRP)],'VariableNames',...
        strsplit('R P'),'RowNames',strsplit('mean std min max')))

end

%% Check parameters
% close all; clc

paracols = [2 13 15:2:21];
paratabs = cell(size(allfits));
for i = 1:length(allfits)    
    paratabs{i} = allfits{i}{:,[2 13 15:2:21]};
    paratabs{i}(:,3:end) = paratabs{i}(:,3:end)./sum(paratabs{i}(:,3:end),2);
%     for j = 1:2
%         tmp = paratabs{i}(:,j);
%         subplot(3,4,(i-1)*4+j)
%         violinplot(tmp);
%         title(sprintf('%.2f±%.2f, range %.2f-%.2f',mean(tmp),std(tmp),min(tmp),max(tmp)))
%     end
%     subplot(3,4,(i-1)*4+[3 4])
%     tmp = paratabs{i}(:,3:end);
%     violinplot(tmp);
end
    
% % Check parameter K
% disp('K')
% tmp = cellfun(@(x) {x(:,1)},paratabs); 
% rmtab = array2table(horzcat(tmp{:}));
% factors = table(categorical(1:3)','VariableName',{'Sess'});
% rm = fitrm(rmtab,sprintf('%s-%s~1',rmtab.Properties.VariableNames{[1 end]}),'WithinDesign',factors);
% [ranovatbl] = ranova(rm, 'WithinModel','Sess'); 
% [H,P,CI,STATS] = ttest(tmp{1},tmp{3});    
% 
% violinplot(horzcat(tmp{:}));hold on
% plot(horzcat(tmp{:})')

% check parameter w
fprintf('\nw_{loc} three-way repeated measures ANOVA\n')
tmp = cellfun(@(x) {x(:,3:6)},paratabs(1:2)); 
t = array2table([tmp{1} tmp{2}]);
factors = table(categorical([2 1 1 2 2 1 1 2])',categorical([1 1 2 2 1 1 2 2])',...
    categorical([1 1 1 1 2 2 2 2])','VariableNames',{'Ecc','VF','Cond'});
rm = fitrm(t,'Var1-Var8~1','WithinDesign',factors);
ranovatbl = ranova(rm,'WithinModel','Ecc*VF*Cond');
etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
disp(ranovatbl(:,[1:5 end]))
mc = multcompare(rm,'Ecc','by','Cond');
disp(mc(mc.Difference>0,:))

% Visualise in violinplots
figure('pos',[-840   439   250*3   159])
for i = 1:2
    subplot(1,3,i)
offs = .4; violinxticks = [1 1 3 3];
v = violinplot(t{:,[2 1 3 4]+(i-1)*4},violinxticks+[-offs offs -offs offs]); %2 1 3 4
for vv = [1 3],v(vv).ViolinColor = ones(1,3)*.2;end
for vv = [2 4],v(vv).ViolinColor = ones(1,3)*.6;end
xticks(unique(violinxticks));xticklabels(strsplit('LVF RVF')); ylabel('w_{loc}')
ylim([0 .8])
box off;
end
subplot(1,3,3)
errorbar(1:2,[mean(reshape(t{:,2:3},np*2,1)) mean(reshape(t{:,[1 4]},np*2,1))],...
    [std(reshape(t{:,2:3},np*2,1)) std(reshape(t{:,[1 4]},np*2,1))],'ko-',...
    'MarkerFaceColor',ones(1,3)*.9)
hold on;
errorbar(1:2,[mean(reshape(t{:,6:7},np*2,1)) mean(reshape(t{:,[5 8]},np*2,1))],...
    [std(reshape(t{:,6:7},np*2,1)) std(reshape(t{:,[5 8]},np*2,1))],'ko-',...
    'MarkerFaceColor',ones(1,3)*.5)
xticks(1:2);xticklabels(strsplit('5\circ 10\circ'));
% lg = legend('Selective','Distributed');
xlim([.5 2.5]); ylim([0 .5]);ylabel('w_{loc}');box off;%xlabel('Eccentricity')

% --- check parameter w in session 2
fprintf('\nS2 w_{loc} two-way repeated measures ANOVA\n')
t = array2table(allfits{3}{:,15:2:21}./sum(allfits{3}{:,15:2:21},2));
factors = table(categorical([2 1 1 2])',categorical([1 1 2 2])',...
    'VariableNames',{'Ecc','VF'});
rm = fitrm(t,'Var1-Var4~1','WithinDesign',factors);
ranovatbl = ranova(rm,'WithinModel','Ecc*VF');
etas = ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
ranovatbl.partialeta = zeros(height(ranovatbl),1);ranovatbl.partialeta(1:2:end) = etas;
disp(ranovatbl(:,[1:5 end]))
% mc = multcompare(rm,'Ecc','by','Cond');
% disp(mc(mc.Difference>0,:))
    
figure('pos',[-440   439   250*1   159])
offs = .4; violinxticks = [1 1 3 3];
plot([0 4],[.25 .25],':','Color',ones(1,3)*.5)
v = violinplot(t{:,[2 1 3 4]},violinxticks+[-offs offs -offs offs]); %2 1 3 4
for vv = [1 3],v(vv).ViolinColor = ones(1,3)*.2;end
for vv = [2 4],v(vv).ViolinColor = ones(1,3)*.6;end
xticks(unique(violinxticks));xticklabels(strsplit('LVF RVF')); ylabel('w_{loc}')
ylim([0 .8])
box off;

fprintf('\n\n***************** End WP2_2_analyse_tvafit *****************\n\n')
    
    
    
    
    