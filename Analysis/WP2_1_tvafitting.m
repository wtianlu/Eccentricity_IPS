% WP2_1_fitTVA
%=======================================================================================
% PROJECT: VF and eccentricity modulations 
% AUTHOR: Lulu Wang
% INSTITUTION: KU Leuven
% AIM: Prepare data for tva fit
% =======================================================================================
% ----------------------------------------------------------------------------------------
% 2020.06.15 Test version with 2 weights instead of 4
function dirs = WP2_1_tvafitting(dirs)


% Load behavioural data
if exist('pinfo','var') == 0, pinfo = readtable([dirs.beh.main '/participant_info.xlsx']);end
p_incl = pinfo.IDCode(~contains(pinfo.Other,'Excluded:')); 
np = length(p_incl);
if exist('bd','var')==0,load([dirs.beh.main '/bd.mat']);end
dirs.tva = '/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_tva1w';% 200615 test
if exist(dirs.tva,'dir')==0,mkdir(dirs.tva);end
nw = 1; % versus 4 weights
disp('Ready!')

%% Prepare data for fitting
clc;
p_refit={};
for p = 1:np
    data = bd.fil_i_em(contains(bd.fil_i_em.pidcode,p_incl{p}),:);
    outputfile = [dirs.tva '/fitS1_em_i_6cond_2w.txt'];
    s = 1;
    fprintf('\n%ss%d_libtva.mat\n',p_incl{p},s)
if exist(sprintf('%ss%d_libtva.mat',p_incl{p},s),'file')>0
    load(sprintf('%ss%d_libtva.mat',p_incl{p},s))
    
else
    % -- Fit session 1
    data1 = data(data.CueI<6,:);
    data1 = data1(data1.session==s,:);
    tvadata = cell(1,height(data1));
    for t = 1:length(tvadata)
        tvadata{t}.condition = str2double(data1.TimI{t})+1;
        tvadata{t}.t = data1.ActTime(t); tvadata{t}.change = data1.Probebool(t) == 2;
        if data1.Probebool(t) == 1, tvadata{t}.probe = data1.Probeloc(t);
        else, tvadata{t}.probe = randi(4);end
    tvadata{t}.response = data1.Response(t)-1; % DK: 2
        
        tvadata{t}.targets = 1:4; tvadata{t}.places = 4; tvadata{t}.display = 1:4;
        tvadata{t}.distractors = []; tvadata{t}.task = 'CD';
    end
    fprintf('Check TVA datafile: %s\nIs CD data: %d\n\n',...
        tvacheckdatafile(tvadata),tvaiscddata(tvadata))
    tvareport(tvadata)
    
    % Construct model
    [theta,tvamodel] = tvainit(tvadata);
    [theta,tvamodel,theta_fix] = tvasculpt(theta,tvamodel,[],'wrep',nw); % N weights: 4 
    [alpha,w,C,svec,v,u0,chdetgm,mu] = tvadeal(tvamodel,1:length(theta));
    theta_fix = nan(u0,1); theta_fix(u0) = 0; theta(u0)=[];
    tvareport(tvadata,tvamodel,theta,theta_fix);
    % Fit model
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.5 100 0.01]); %shave K
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.1  1000 0.01]);
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.02 1000 0.01]);
    % Save data
    tvareport(tvadata,tvamodel,theta,theta_fix);
    if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta,theta_fix);end
    tvalpr(outputfile,sprintf('%ss%d',p_incl{p},s),tvadata,tvamodel,theta,theta_fix);
    save(sprintf('%ss%d_libtva.mat',p_incl{p},s),'tvadata','tvamodel','theta','theta_fix');
    fprintf('Finished fitting %ss%d\n\n',p_incl{p},s)    
end

% --- Check fit performance -> fix C?
[Rsq,oo,pp] = tvarsq(tvadata,tvamodel,theta,theta_fix);[R,P] = corrcoef(oo,pp);
[humtheta,sigma,ll,pM, pH,normg,H,g,U,rcondU] = tvadiag(theta,tvamodel,tvadata,theta_fix);
fprintf('pH = %d Rsq: %.3f R: %.3f p: %.3f\n',pH,Rsq,R(1,2),P(1,2))
% pause;
if or(P(1,2)>0.049,pH>0)
    p_refit = [p_refit;p_incl(p)];
    continue;
end
theta1 = theta;

% -- Fit session 2
s = 2; fprintf('%ss%d_libtva.mat\n',p_incl{p},s);
if exist(sprintf('%ss%d_libtva.mat',p_incl{p},s),'file')==0
% for p = 1
%     load(sprintf('%ss%d_libtva.mat',p_incl{p},s)) % previous session
    
    outputfile = [dirs.tva '/fitS2_em_i_1cond.txt'];

    data1 = data(data.session==s,:);
    tvadata = cell(1,height(data1));
    for t = 1:length(tvadata)
        tvadata{t}.t = data1.ActTime(t); tvadata{t}.change = data1.Probebool(t) == 2;
        if data1.Probebool(t) == 1, tvadata{t}.probe = data1.Probeloc(t);
        else, tvadata{t}.probe = randi(4);end
        tvadata{t}.response = data1.Response(t)-1; % DK: 2

        tvadata{t}.condition = 1;
        tvadata{t}.targets = 1:4; tvadata{t}.places = 4; tvadata{t}.display = 1:4;
        tvadata{t}.distractors = []; tvadata{t}.task = 'CD';
    end
    fprintf('Check TVA datafile: %s\nIs CD data: %d\n\n',...
        tvacheckdatafile(tvadata),tvaiscddata(tvadata))
    tvareport(tvadata)
    % Construct model
    [theta,tvamodel] = tvainit(tvadata);
    [theta,tvamodel,theta_fix] = tvasculpt(theta,tvamodel,[],'wrep',nw); % N weights: 4
    [alpha,w,C,svec,v,u0,chdetgm,mu] = tvadeal(tvamodel,1:length(theta));
    theta_fix = nan(u0,1); theta_fix(u0) = 0; theta_fix(C) = theta1(C);
    theta(u0)=[]; theta(C) = [];
    tvareport(tvadata,tvamodel,theta,theta_fix);
    % Fit model
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.5 100 0.01]); %shave K
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.1  1000 0.01]);
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.02 1000 0.01]);
    % Save data
    tvareport(tvadata,tvamodel,theta,theta_fix);
    if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta,theta_fix);end
    tvalpr(outputfile,sprintf('%ss%d',p_incl{p},s),tvadata,tvamodel,theta,theta_fix);
    save(sprintf('%ss%d_libtva.mat',p_incl{p},s),'tvadata','tvamodel','theta','theta_fix');
    fprintf('Finished fitting %ss%d\n\n',p_incl{p},s)   
else
    load(sprintf('%ss%d_libtva.mat',p_incl{p},s))
end
% --- Check fit performance 
[humtheta,sigma,ll,pM, pH,normg,H,g,U,rcondU] = tvadiag(theta,tvamodel,tvadata,theta_fix);
fprintf('pH = %d\n',pH)
% pause;

% -- Fit baseline
s = 1; fprintf('%s_bl_s%d_libtva.mat\n',p_incl{p},s)
if exist(sprintf('%s_bl_s%d_libtva.mat',p_incl{p},s),'file')==0
    outputfile = [dirs.tva '/fitS1baseline_em_i_6cond.txt'];
    data1 = data(data.session==s,:); data1 = data1(data1.CueI==6,:);
    tvadata = cell(1,height(data1));
    for t = 1:length(tvadata)
        tvadata{t}.t = data1.ActTime(t); tvadata{t}.change = data1.Probebool(t) == 2;
        if data1.Probebool(t) == 1, tvadata{t}.probe = data1.Probeloc(t);
        else, tvadata{t}.probe = randi(4);end
        tvadata{t}.response = data1.Response(t)-1;

        tvadata{t}.condition = str2double(data1.TimI{t})+1;
        tvadata{t}.targets = 1:4; tvadata{t}.places = 4; tvadata{t}.display = 1:4;
        tvadata{t}.distractors = []; tvadata{t}.task = 'CD';
    end
    fprintf('Check TVA datafile: %s\nIs CD data: %d\n\n',...
        tvacheckdatafile(tvadata),tvaiscddata(tvadata))
    tvareport(tvadata)
    % Construct model
    [theta,tvamodel] = tvainit(tvadata);
    [theta,tvamodel,theta_fix] = tvasculpt(theta,tvamodel,[],'wrep',nw); % N weights: 4
    [alpha,w,C,svec,v,u0,chdetgm,mu] = tvadeal(tvamodel,1:length(theta));
    theta_fix = nan(u0,1); theta_fix(u0) = 0; theta_fix(C) = theta1(C);
    theta(u0)=[]; theta(C) = [];
    tvareport(tvadata,tvamodel,theta,theta_fix);
    % Fit model
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.5 100 0.01]); %shave K
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.1  1000 0.01]);
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.02 1000 0.01]);
    % Save data
    tvareport(tvadata,tvamodel,theta,theta_fix);
    if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta,theta_fix);end
    tvalpr(outputfile,sprintf('%s_bl_s%d',p_incl{p},s),tvadata,tvamodel,theta,theta_fix);
    save(sprintf('%s_bl_s%d_libtva.mat',p_incl{p},s),'tvadata','tvamodel','theta','theta_fix');
    fprintf('Finished fitting %ss%d_bl\n\n',p_incl{p},s)  
else
    load(sprintf('%s_bl_s%d_libtva.mat',p_incl{p},s))
end
[Rsq,oo,pp] = tvarsq(tvadata,tvamodel,theta,theta_fix);[R,P] = corrcoef(oo,pp);
[humtheta,sigma,ll,pM, pH,normg,H,g,U,rcondU] = tvadiag(theta,tvamodel,tvadata,theta_fix);
fprintf('pH = %d Rsq: %.3f R: %.3f p: %.3f\n\n',pH,Rsq,R(1,2),P(1,2))

end


% No good fits

% p_refit = strsplit('TE07');
s = 1;
for p = 1:length(p_refit)
fprintf('\n%ss%d_libtva_fixC.mat\n',p_refit{p},s)
% -- re-fit S1    
outputfile = [dirs.tva '/fitS1_em_i_6cond_fixC.txt'];
if exist(sprintf('%ss%d_libtva_fixC.mat',p_refit{p},s),'file')==0
    load(sprintf('%ss%d_libtva.mat',p_refit{p},s))
    clear theta theta_fix tvamodel; % keep tvadata
    % Construct model
    [theta,tvamodel] = tvainit(tvadata);
    [theta,tvamodel,theta_fix] = tvasculpt(theta,tvamodel,[],'wrep',nw); % N weights: 4
    [alpha,w,C,svec,v,u0,chdetgm,mu] = tvadeal(tvamodel,1:length(theta));
    theta(u0)=[]; theta(C) = [];
    theta_fix = nan(u0,1); theta_fix(u0) = 0; theta_fix(C) = -2.52573;%tvafixer('+',60);
    tvareport(tvadata,tvamodel,theta,theta_fix); % 50: -2.9957; 60: -2.8134; 80: -2.52573
    % Fit model
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.5 100 0.01]); %shave K
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.1  1000 0.01]);
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.02 1000 0.01]);
    % Save data
    tvareport(tvadata,tvamodel,theta,theta_fix);
    if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta,theta_fix);end
    tvalpr(outputfile,sprintf('%ss%d',p_refit{p},s),tvadata,tvamodel,theta,theta_fix);
    save(sprintf('%ss%d_libtva_fixC.mat',p_refit{p},s),'tvadata','tvamodel','theta','theta_fix');
    fprintf('Finished fitting %ss%d\n\n',p_refit{p},s)   
else
    load(sprintf('%ss%d_libtva_fixC.mat',p_refit{p},s))
end
[Rsq,oo,pp] = tvarsq(tvadata,tvamodel,theta,theta_fix);[R,P] = corrcoef(oo,pp);
fprintf('Rsq: %.3f R: %.3f p: %.3f\n',Rsq,R(1,2),P(1,2))
pause;
end

% Check s1 fit without fixing t0
fits1data = dir([dirs.main '/Data_tva/s1_fixt0/TE*s1_libtva.mat']);
outputfile = [dirs.main '/Data_tva/fitS1_em_i_6cond_fixnone.txt']; s=1;
for p = 2:length(fits1data)
    load([fits1data(p).folder '/' fits1data(p).name],'tvadata')
    [theta,tvamodel] = tvafit(tvadata,4);
    if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta);end
    tvalpr(outputfile,sprintf('%ss%d',p_incl{p},s),tvadata,tvamodel,theta);
    save(sprintf('%s/Data_tva/s1_fixnone/%ss%d_libtva_fixnone.mat',dirs.main,p_incl{p},s),...
        'tvamodel','theta');
    fprintf('Finished fitting %ss%d\n\n',p_incl{p},s)       
end

% Fit session 1 baseline with final S1 fit parameters
clc;
dirs.tvafit1 = [dirs.main '/Data_tva/s1_final'];
for p = 1:length(p_incl)
    
    % Load t0 and C from s1
    fit1mat = dir([dirs.tvafit1 '/' p_incl{p} 's1_*.mat']);
    if length(fit1mat)~=1,disp(fit1mat);pause;end
    disp([fit1mat.folder '/' fit1mat.name])
    
    fitout = sprintf('%s/Data_tva/s1bl_fixt0C/%s_bl_s%d_libtva.mat',dirs.main,p_incl{p},s);
    if exist(fitout,'file')==0
    load([fit1mat.folder '/' fit1mat.name])
    theta1 = theta; theta_fix1 = theta_fix; theta1 = tvathetacombine(theta1,theta_fix1);
    clear tvadata tvamodel theta theta_fix;
    
    % -- Fit baseline data
    s = 1; outputfile = [dirs.tva '/fitS1baseline_em_i_6cond.txt'];
    fitblmat = dir([dirs.main '/Data_tva/s1bl_tmp/' p_incl{p} '_bl_s1_*.mat']);
    if length(fitblmat)~=1,disp(fitblmat);pause;end
    load([fitblmat.folder '/' fitblmat.name],'tvadata');
    % Construct model
    [theta,tvamodel] = tvainit(tvadata);
    [theta,tvamodel,theta_fix] = tvasculpt(theta,tvamodel,[],'wrep',nw); % N weights: 4
    [~,w,C,~,~,u0,~,~] = tvadeal(tvamodel,1:length(theta));
    theta_fix = nan(u0,1); theta_fix(u0) = 0; theta_fix(C) = theta1(C);
    theta(u0)=[]; theta(C) = [];
    % Fit model
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.5 100 0.01]); %shave K
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.1  1000 0.01]);
    [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.02 1000 0.01]);
    % Save data
    tvareport(tvadata,tvamodel,theta,theta_fix);
    if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta,theta_fix);end
    tvalpr(outputfile,sprintf('%s_bl_s%d',p_incl{p},s),tvadata,tvamodel,theta,theta_fix);
    save(fitout,'tvadata','tvamodel','theta','theta_fix');
    fprintf('\nFinished fitting %ss%d_bl\n\n',p_incl{p},s)  
    else
        load(fitout)
    end
    
    
%     tvareport(tvadata,tvamodel,theta,theta_fix);
[Rsq,oo,pp] = tvarsq(tvadata,tvamodel,theta,theta_fix);[R,P] = corrcoef(oo,pp);
[humtheta,sigma,ll,pM, pH,normg,H,g,U,rcondU] = tvadiag(theta,tvamodel,tvadata,theta_fix);
fprintf('\npH = %d Rsq: %.3f R: %.3f p: %.3f\n\n',pH,Rsq,R(1,2),P(1,2))

end

% Fit session 2 with final S1 fit parameters
clc;
outputfile = [dirs.tva '/fitS2_em_i_6cond_fixt0C.txt'];
for p = 1:np
    % -- Fit session 2 data
    s = 2;
    fitout = sprintf('%s/Data_tva/s2_fixt0C/%ss%d_libtva.mat',dirs.main,p_incl{p},s); 
    disp(fitout)
    if exist(fitout,'file')~=0
        load(fitout)
    else
        % load session 1 model
        fit1mat = dir([dirs.tvafit1 '/' p_incl{p} 's1_*.mat']);
        if length(fit1mat)~=1,disp(fit1mat);pause;end
        load([fit1mat.folder '/' fit1mat.name])
        theta1 = tvathetacombine(theta,theta_fix);
        clear theta theta_fix tvamodel tvadata
        
        % load session 2 tvadata
        fit2mat = dir([dirs.main '/Data_tva/s2_fixt0/' p_incl{p} 's2_*.mat']);
        if length(fit2mat)~=1,disp(fit2mat);pause;end
        load([fit2mat.folder '/' fit2mat.name],'tvadata');
        
        % Construct model
        [theta,tvamodel] = tvainit(tvadata);
        [theta,tvamodel,theta_fix] = tvasculpt(theta,tvamodel,[],'wrep',nw); % N weights: 4
        [~,w,C,~,~,u0,~,~] = tvadeal(tvamodel,1:length(theta));
        theta_fix = nan(u0,1); theta_fix(u0) = 0; theta_fix(C) = theta1(C);
        theta(u0)=[]; theta(C) = [];
        % Fit model
        [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.5 100 0.01]); 
        [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.1  1000 0.01]);
        [theta,tvamodel,change] = tvashave(theta,tvamodel,tvadata,theta_fix,[0.02 1000 0.01]);
        % Save data
        tvareport(tvadata,tvamodel,theta,theta_fix);
        if exist(outputfile,'file')==0,tvalpr(outputfile,'',tvadata,tvamodel,theta,theta_fix);end
        tvalpr(outputfile,sprintf('%s_bl_s%d',p_incl{p},s),tvadata,tvamodel,theta,theta_fix);
        save(fitout,'tvadata','tvamodel','theta','theta_fix');
        fprintf('\nFinished fitting %ss%d_bl\n\n',p_incl{p},s)  
    
    end
    
    [humtheta,sigma,ll,pM, pH,normg,H,g,U,rcondU] = tvadiag(theta,tvamodel,tvadata,theta_fix);
    fprintf('\npH = %d\n\n',pH)
%     pause;
end

