%=======================================================================================
% WP2_1_prepdata_beh
%=======================================================================================
% PROJECT:      WP2 - VF and eccentricity modulations 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Prepare all behavioiral data across all participants for pre-processing
% OUTPUT:       Beh data organisedin /Data_proc_beh
% -------------------------------------------------------------------------
% 2020.05.30 Half
% -------------------------------------------------------------------------

function dirs = WP2_1_prepdata_beh(dirs,qc)
% Initialisation
if nargin<2, qc=false; end
fprintf('\n\n***************** Start WP2_1_prepdata_beh *****************\n\n')

%% Read participant information
pinfo = readtable([dirs.beh.main '/participant_info.xlsx']);
p_incl = pinfo.IDCode(~contains(pinfo.Other,'Excluded:')); 
np = length(p_incl);
fprintf('Included participants: %d\n',np);fprintf('%s\t',p_incl{:});fprintf('\n')

%% Get raw behavioural data

sesfn = strsplit('BT fT'); validitycode = {1,2,3,4,[1 4],[2 3]};

dirs.beh.mat = [dirs.beh.main '/EccVF_rawdata.mat'];
if exist(dirs.beh.mat,'file')==0
    % Get filepaths
    dirs.raw.beh = [dir([dirs.raw.main '/TE*/rawdata/P*_BTask.txt']);
        dir([dirs.raw.main '/TE*/rawdata/P*_fTask.txt'])];
    dirs.raw.beh = fullfile({dirs.raw.beh.folder}',{dirs.raw.beh.name}');
    dirs.raw.beh = dirs.raw.beh(contains(dirs.raw.beh,p_incl));
    % Read raw behavioural data
    rds = pinfo(~contains(pinfo.Other,'Excluded:'),[1:4,6]);
    for p = 1:length(p_incl)
        for s = 1:2
        % load session data
            fn{s} = dirs.raw.beh{and(contains(dirs.raw.beh,p_incl{p}),...
                contains(dirs.raw.beh,sesfn{s}))};
            disp(['Loading ' fn{s}]);
            fid=fopen(fn{s}); tline = fgetl(fid); data_raw = cell(0,1);
            while ischar(tline),data_raw{end+1,1} = tline;tline = fgetl(fid);end
            fclose(fid); clear tline;
        % process data
            trialidx = cellfun(@(x) ~isempty(strfind(x,'T:')),data_raw);
            fprintf('%d trials found in %s\n',sum(trialidx),fn{s})
            data_tmp = cellfun(@(x) {strsplit(x,'\t')},data_raw(trialidx));
            % separate trials
            data_trials = num2cell(nan(size(data_tmp,1),max(cellfun(@length,data_tmp))));
            for l=1:size(data_tmp,1), data_trials(l,1:length(data_tmp{l}))=data_tmp{l}; end
            data_trials(:,[1:6 9 10]) = num2cell(cellfun(@str2double,data_trials(:,[1:6 9 10])));
            datatbl = cell2table(data_trials,'VariableNames',...
                strsplit('TrialI CueI PrbI TimI Frames Config Targets Probe Onset ActTime RawResp RT Trialtimings'));
            % clean data
            idx = cellfun(@length,datatbl.RawResp)>2;
            Response = nan(height(datatbl),1);
            Response(idx) = cellfun(@(x) str2double(x(3)),datatbl.RawResp(idx));
            Validity = nan(height(datatbl),1); Probebool = nan(height(datatbl),1);
            for i = 1:length(Validity)
                if contains(datatbl.Targets{i}(3:end),datatbl.Probe{i}(3:end))==0
                    Validity(i)=2; %NP
                    Probebool(i) = 2;
                else
                    Probebool(i) = 1;
                    if datatbl.CueI(i)==6,Validity(i)=-1;%DISTR
                    elseif ismember(find(datatbl.Targets{i}(3:end)==datatbl.Probe{i}(3:end)),...
                            validitycode{1+datatbl.CueI(i)}), Validity(i)=1;%VALID
                    else, Validity(i) = 0; %INVALID
                    end
                end
            end
            Correct = nan(height(datatbl),1);
        end      
        % save data
        rds.fn = fn;
        
    end
    save(dirs.beh.mat,'rds')
else
    load(dirs.beh.mat)
end

disp([dirs.beh '/EccVF_rawdata.txt created!'])









