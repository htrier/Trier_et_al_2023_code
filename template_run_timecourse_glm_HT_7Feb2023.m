%% 1. Select ROIs and model
clear all; fclose all; clc;
warning('on','all');
GLM = 3; % 1=basic or 2=PPI or 3=PPI combining forage and check actions from vigilance phase
stat = 1;
textout = 0;
plot_session = 1;
plot_session_group = 0;
upsample = 20;
pre_win = 2; % 5
drive_name = 'Seagate_storage';

epi_dir = strcat(['/Volumes/',drive_name,'/Hailey_data_preproc/fMRI_data/preprocess/']);
behav_dir = strcat(['/Volumes/',drive_name,'/Hailey_data_preproc/fMRI_data/']);
timecourse_dir = strcat(['/Volumes/',drive_name,'/Hailey_data_preproc/fMRI_data/group_analyses/masks/selected_masks_for_ROI_analysis/2mm_masks_22Feb2023/masks/']);
% timecourse_dir = strcat(['/Volumes/',drive_name,'/Hailey_data_preproc/fMRI_data/group_analyses/masks/selected_masks_for_ROI_analysis/1mm/masks/']);
feat_name = 'L1stats_27Apr2022_noTempDerivs_motionRegressors_NoInterestContrasts';

animals = {'sub103','sub105','sub109','sub111','sub117','sub119','sub121','sub123','sub125','sub127',...
    'sub129','sub131','sub204','sub208','sub212','sub214','sub216','sub218','sub222','sub224','sub228','sub230','sub232'};

% Relevant ROIs (27 Feb 2023): 'ACC_cluster_cope7','DRN_sub_PAG_edit_sub_SC_edit_ero','HB','insular_cortex_cope7',...
%     'motor_sphere_rad2','NAc','SN_Pauli'...
%     'VTA_Nima_Pauli_sub_SN_Pauli','PAG_ero','SC_sub_PAG','striatal_sphere_rad4_dil'

% Threat level + check analyses (2s pre-window)
% roi= {'HB','DRN_sub_PAG_edit_sub_SC_edit_ero','HB','DRN_sub_PAG_edit_sub_SC_edit_ero',...
%     'VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero','SN_Pauli'};
% seed_roi = {'ACC_cluster_cope7','ACC_cluster_cope7','insular_cortex_cope7','insular_cortex_cope7','HB','HB','HB'};

% Threat level + check analyses (2s pre-window)
% roi = {'DRN_sub_PAG_edit_sub_SC_edit_ero','HB'};
% seed_roi = {'insular_cortex_cope7','insular_cortex_cope7'};
% tc_action = 'forages'; % if focusing on only one action type, indicate here (for running a model with no action regressor)

% Threat level + check analyses (5s pre-window)
% roi = {'HB'};
% seed_roi = {'ACC_cluster_cope7'};
% tc_action = 'checks'; % if focusing on only one action type, indicate here (for running a model with no action regressor)

% % Threat level + Successful/unsuccessful check analyses
% roi = {'SN_Pauli','VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero'};
% seed_roi = {'HB','HB','HB'};

% Reward + forage analyses (2s pre-window)
% roi = {'striatal_sphere_rad4_dil','SN_Pauli',...
%     'VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero',...
%     'striatal_sphere_rad4_dil','DRN_sub_PAG_edit_sub_SC_edit_ero','SN_Pauli','VTA_Nima_Pauli_sub_SN_Pauli'};
% seed_roi = {'precentral_constant_rad2_7Feb','precentral_constant_rad2_7Feb',...
%     'precentral_constant_rad2_7Feb','precentral_constant_rad2_7Feb'...
%   'ACC_cluster_cope7','ACC_cluster_cope7','ACC_cluster_cope7','ACC_cluster_cope7'};
roi={'DRN_sub_PAG_edit_sub_SC_edit_ero'};
seed_roi = {'precentral_constant_rad2_7Feb'};
% tc_action='forages';

% Reward + forage analyses (2s pre-window): One action type only
% roi = {'DRN_sub_PAG_edit_sub_SC_edit_ero'};
% seed_roi = {'ACC_cluster_cope7'};
% tc_action = 'forages';

% Threat level + forage (5s pre-window)
% roi = {'HB'};
% seed_roi = {'ACC_cluster_cope7'};

% Threat level + forage analyses (2s pre-window)
% roi = {'SN_Pauli','VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero'};
% seed_roi = {'HB','HB','HB'};

% Reward + check analyses (2s pre-window)
% roi = {'SN_Pauli','VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero'...
%     'striatal_sphere_rad4_dil','DRN_sub_PAG_edit_sub_SC_edit_ero',...
%     'striatal_sphere_rad4_dil','DRN_sub_PAG_edit_sub_SC_edit_ero'};
% seed_roi = {'HB','HB','HB',...
%     'precentral_constant_rad2_7Feb','precentral_constant_rad2_7Feb',...
%     'ACC_cluster_cope7','ACC_cluster_cope7'};

% Analyses to compare SN/VTA interactions with results shown in 3a and 3b
% roi = {'SN_Pauli','VTA_Nima_Pauli_sub_SN_Pauli','SN_Pauli','VTA_Nima_Pauli_sub_SN_Pauli'};
% seed_roi = {'ACC_cluster_cope7','ACC_cluster_cope7','insular_cortex_cope7','insular_cortex_cope7'};

%BehavAll = readtable('/Volumes/MyPassport/Hailey_data_preproc/fMRI_data/all_subjects_behavior.csv');%[behav_dir,'all_subjects_behavior.csv']);
if pre_win == 2
    load(strcat('/Volumes/',drive_name,'/Hailey_data_preproc/fMRI_data/all_subjects_behav_updatedEpochedData_upsample20_minusHalfTR_15May2023.mat'));
%     load(strcat('/Volumes/MyPassport/Hailey_data_preproc/fMRI_data/all_subjects_behav_updatedEpochedData_upsample',num2str(upsample),'_minusHalfTR_7Feb2023_RewardLags.mat'));
elseif pre_win == 5
    load(strcat('/Volumes/',drive_name,'/Hailey_data_preproc/fMRI_data/all_subjects_behav_updatedEpochedData_upsample20_minusHalfTR_5sPreWindow.mat'));
end

% Choose model
% Pre-discovery phase
% currModel = {'switchToCheck','timepressure','time','reward','seedTC','constant','PPI1_seedTCXswitchToCheck','PPI2_seedTCXtimepressure','PPI3_switchToCheckXtimepressure','PPI4_seedTCXswitchToCheckXtimepressure'};
currModel = {'switchToForage','reward','time','timepressure','seedTC','constant','PPI1_seedTCXswitchToForage','PPI2_seedTCXreward','PPI3_switchToForageXreward','PPI4_seedTCXswitchToForageXreward'};
% currModel = {'switchToForage','rewardAvg8sLag','time','timepressure','seedTC','constant','PPI1_seedTCXswitchToForage','PPI2_seedTCXrewardAvg8sLag','PPI3_switchToForageXrewardAvg8sLag','PPI4_seedTCXswitchToForageXrewardAvg8sLag'};
% currModel = {'switchToForage','rewardAvg4sLag','time','timepressure','seedTC','constant','PPI1_seedTCXswitchToForage','PPI2_seedTCXrewardAvg4sLag','PPI3_switchToForageXrewardAvg4sLag','PPI4_seedTCXswitchToForageXrewardAvg4sLag'};
% currModel = {'switchToForage','timepressure','time','reward','seedTC','constant','PPI1_seedTCXswitchToForage','PPI2_seedTCXtimepressure','PPI3_switchToForageXtimepressure','PPI4_seedTCXswitchToForageXtimepressure'};
% currModel = {'switchToCheck','reward','time','timepressure','seedTC','constant','PPI1_seedTCXswitchToCheck','PPI2_seedTCXreward','PPI3_switchToCheckXreward','PPI4_seedTCXswitchToCheckXreward'};
% currModel = {'timepressure','time','reward','seedTC','constant','PPI2_seedTCXtimepressure'};
% currModel = {'reward','time','timepressure','seedTC','constant','PPI2_seedTCXreward'};

% Post-discovery phase
% a = {'switchToForage','reward','time','seedTC','constant','PPI1_seedTCXswitchToForage','PPI2_seedTCXreward','PPI3_switchToForageXreward','PPI4_seedTCXswitchToForageXreward'};
% currModel = {'switchToCheck','proximity','time','reward','seedTC','constant','PPI1_seedTCXswitchToCheck','PPI2_seedTCXproximity','PPI3_switchToCheckXproximity','PPI4_seedTCXswitchToCheckXproximity'};

% Map regressor names to nice labels for final figures
regnames = {'switchToCheck','timepressure','time','reward','seedTC','constant',...
    'PPI1_seedTCXswitchToCheck','PPI2_seedTCXtimepressure',...
    'PPI3_switchToCheckXtimepressure','PPI4_seedTCXswitchToCheckXtimepressure',...
    'switchToForage','PPI1_seedTCXswitchToForage','PPI2_seedTCXreward',...
    'PPI3_switchToForageXreward','PPI4_seedTCXswitchToForageXreward'...
    'PPI3_switchToForageXtimepressure','PPI4_seedTCXswitchToForageXtimepressure',...
    'PPI3_switchToCheckXreward','PPI4_seedTCXswitchToCheckXreward',...
    'PPI3_switchToCheckXproximity','PPI2_seedTCXproximity','PPI4_seedTCXswitchToCheckXproximity',...
    'proximity'};
reglabels = {'Check','Threat level','Time','Reward','seedTC','Constant',...
    'seedTC*Check','seedTC*Threat level','Check*Threat level','seedTC*Check*Threat level',...
    'Forage','seedTC*Forage','seedTC*Reward','Forage*Reward','seedTC*Forage*Reward',...
    'Forage*Threat level','seedTC*Forage*Threat level','Check*Reward','seedTC*Check*Reward',...
    'Check*Proximity','seedTC*Proximity','seedTC*Check*Proximity','Proximity'};
regLabelMap = containers.Map(regnames, reglabels);
clear regnames, reglabels;

% Map ROI names to nice labels
roinames = {'striatal_sphere_rad4_dil','SN_Pauli',...
    'VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero',...
    'precentral_constant_rad2_7Feb','ACC_cluster_cope7','insular_cortex_cope7','HB'};
roilabels = {'Striatum','SN','VTA','DRN','Precentral gyrus','ACC','AI','LHb'};
roiLabelMap = containers.Map(roinames, roilabels);
clear roinames, roilabels;

%% 2. Run GLM
% Possible conditions: vigilance_firstChecks_constant /
% vigilance_successfulChecks_outcome / vigilance_unsuccessfulChecks /
% monitoring_firstChecks_constant (to get both first checks & forages from
% post-discovery phase)
cond = 'vigilance_firstChecks_constant';
% Get data for this condition - vigilance_successfulChecks_outcome is a special case
if contains(cond,'vigilance_successfulChecks_outcome')
    firstChecksData = updatedBehav(strcmp(updatedBehav.action, 'allChecks') & strcmp(updatedBehav.phase, 'vigilance') & updatedBehav.first_check==1 & updatedBehav.successful_check==1, :);
    firstForagesData = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'vigilance'), :);

elseif contains(cond,'vigilance_unsuccessfulChecks')
    firstChecksData = updatedBehav(strcmp(updatedBehav.action, 'allChecks') & strcmp(updatedBehav.phase, 'vigilance') & updatedBehav.first_check==1 & updatedBehav.successful_check==0, :);
    firstForagesData = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'vigilance'), :);

elseif contains(cond,'vigilance_firstChecks_constant')
    firstChecksData = updatedBehav(strcmp(updatedBehav.action, 'allChecks') & strcmp(updatedBehav.phase, 'vigilance') & updatedBehav.first_check==1, :);
    firstForagesData = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'vigilance'), :);

elseif contains(cond,'monitoring_firstChecks_constant')
    firstChecksData = updatedBehav(strcmp(updatedBehav.action, 'firstChecks') & strcmp(updatedBehav.phase, 'monitoring'), :);
    firstForagesData = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'monitoring'),:);
end

for r = 1:length(roi)  % r=1;
    z=0;
    for a = 1:length(animals)
        animal = animals{a};
        z = z+1;
    
        subj_firstChecksData = firstChecksData(strcmp(firstChecksData.subjectNumber, animal), :);
        subj_firstForagesData = firstForagesData(strcmp(firstForagesData.subjectNumber, animal), :);
    
        % Run time course
        %main
        if contains(currModel(1),'switch')
            REG.switchToCheck = [repmat(1,[size(subj_firstChecksData,1),1]); repmat(-1,[size(subj_firstForagesData,1),1])]; %1=check, -1=forage
            REG.switchToForage = [repmat(-1,[size(subj_firstChecksData,1),1]); repmat(1,[size(subj_firstForagesData,1),1])]; %-1=check, 1=forage
            REG.constant = ones (length (REG.switchToCheck),1);
            REG.time = [subj_firstChecksData.('onset'); subj_firstForagesData.('onset')]; % checks, then forages
            REG.reward = [subj_firstChecksData.('reward'); subj_firstForagesData.('reward')]; % checks, then forages
            REG.timepressure = [subj_firstChecksData.('timePressure'); subj_firstForagesData.('timePressure')];
            REG.posUncertainty = [subj_firstChecksData.('posUncertainty'); subj_firstForagesData.('posUncertainty')];
            REG.proximity = [subj_firstChecksData.('proximity'); subj_firstForagesData.('proximity')];
        else
            % Do not combine checks and forages
            if strcmp(tc_action,'checks')
                REG.constant = repmat(1,[size(subj_firstChecksData,1),1]);
                REG.time = [subj_firstChecksData.('onset')];
                REG.reward = [subj_firstChecksData.('reward')];
                REG.timepressure = [subj_firstChecksData.('timePressure')];
                REG.posUncertainty = [subj_firstChecksData.('posUncertainty')];
            elseif strcmp(tc_action,'forages')
                REG.constant = repmat(1,[size(subj_firstForagesData,1),1]);
                REG.time = [subj_firstForagesData.('onset')]; % checks, then forages
                REG.reward = [subj_firstForagesData.('reward')]; % checks, then forages
                REG.timepressure = [subj_firstForagesData.('timePressure')];
                REG.posUncertainty = [subj_firstForagesData.('posUncertainty')];
            end
        end

        % load time-series data
        if pre_win==5
            data_identifier = '_minusHalfTR_5sPreWindow_epoched';
        else
            data_identifier = '_minusHalfTR_epoched';
        end
    
        if strcmp(cond, 'monitoring_firstChecks_constant')
            seed_TC_checks = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',seed_roi{r},'_monitoring_firstChecks_constant_upsample',num2str(upsample),data_identifier])]);
            seed_TC_forages = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',seed_roi{r},'_monitoring_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
    
            roi_TC_checks  = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',roi{r},'_monitoring_firstChecks_constant_upsample',num2str(upsample),data_identifier])]);
            roi_TC_forages  = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',roi{r},'_monitoring_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
        else
            seed_TC_checks = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',seed_roi{r},'_',cond,'_upsample',num2str(upsample),data_identifier])]);
            seed_TC_forages = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',seed_roi{r},'_vigilance_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
            
            roi_TC_checks  = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',roi{r},'_',cond,'_upsample',num2str(upsample),data_identifier])]);
            roi_TC_forages  = load([timecourse_dir,animal,'/',feat_name,'/',(['tc_',roi{r},'_vigilance_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
        end
        
        if contains(currModel(1),'switch')
            seed_TC = [seed_TC_checks.trial_data; seed_TC_forages.trial_data];
            roi_TC = [roi_TC_checks.trial_data; roi_TC_forages.trial_data];
        else
            % Do not combine checks and forages
            if strcmp(tc_action,'checks')
                seed_TC = [seed_TC_checks.trial_data];
                roi_TC = [roi_TC_checks.trial_data];
            elseif strcmp(tc_action,'forages')
                seed_TC = [seed_TC_forages.trial_data];
                roi_TC = [roi_TC_forages.trial_data];
            end

        end
        clear seed_TC_checks seed_TC_forages roi_TC_checks roi_TC_forages
    
        % run GLM
        betaOut = NaN(size(currModel,2),size(roi_TC,2));
        for i = 1:size(roi_TC,2)
    
            % load ROI time-series data               
            REG.seedTC  = seed_TC(:,i);
    
            % Make regressor design matrix (don't forget constant!)
            dmat = [];
            for c = 1:length(currModel)
                if ~contains(currModel{c},'PPI')
%                     disp(strcat(['size dmat: ',num2str(size(dmat,1)),' x ',num2str(size(dmat,2)),' / size REG (',currModel{c},'): ',num2str(size(REG.(currModel{c}),1)),' x ',num2str(size(REG.(currModel{c}),2))]));
                    dmat = [dmat REG.(currModel{c})];
                end
            end
    
            % remove trials with no response
            dmat(isnan(dmat(:,2)),:)=[];
    
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1))); % normalize binary action var
%             dmat(:,2:(size(dmat,2)-1)) = zscore(dmat(:,2:(size(dmat,2)-1))); % leave action var unnormalized
            
            % Add PPI regressors
            for ppi = 1:length(currModel)
                if contains(currModel{ppi},'PPI')
                    % Determine vars that go into this PPI interaction
                    ppivars = split(currModel{ppi},'_');
                    ppivars = split(ppivars{2},'X');
                    if length(ppivars)==2
                        REG.(currModel{ppi}) = zscore(REG.(ppivars{1}) .* REG.(ppivars{2}));
                    elseif length(ppivars)==3
                        REG.(currModel{ppi}) = zscore(REG.(ppivars{1}) .* REG.(ppivars{2}) .* REG.(ppivars{3}));
                    end
                    dmat = [dmat, REG.(currModel{ppi})];
                end
            end
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
    
            % beta X time output
            betaOut(:,i) = ols(roi_TC(:,i),dmat,contrasts);
            clear dmat contrasts REG.checkNotForage REG.seedTC REG.time REG.constant REG.PPI1_seedTCXcheckNotForage
        end
    
        allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',strcat([cond,'_plus_firstForages'])))(:,:,z) = betaOut;
        clear betaOut
    end
end
%% 3. Test, plot, and save stats
tableVarNames = {'Condition','Regressor','ROI','Seed','Uncorrected_pvalue','Significant',...
    'T_stats','df','Avg_peak_sec','Peak_search_window','Mean','SD','p_wilcoxon','zval_wilcoxon','signedrank_wilcoxon','Full_model',...
    'Peak_vals'}; % 'HolmBonf_corrected_pvalue'
allstats = cell2table(cell(0,17), 'VariableNames', tableVarNames);
format shortG

clearvars window LOOT peak
numsession = 1:numel(animals);
post_win = 15;
TR = 1.962;
nsamples = round(((pre_win+post_win)./TR)*upsample);
time = -pre_win:(pre_win+post_win)/nsamples:post_win;
makeplot = 1;
if makeplot==1
    clr = [0.8,0,0];
    plt_colour = {clr};
end

if pre_win==2
    start_inds = [22,72]; % early [0.0636s, 4.9769s], late [4.9769s, 9.9884s]
    end_inds = [72,123]; % early [0.0636s, 4.9769s], late [4.9769s, 9.9884s]
elseif pre_win==5
    start_inds = [1,103]; % if 5s prewindow: [1,103] for -5-5s 
    end_inds = [103,154]; % if 5s prewindow: [103,154] for 5-10s
end

for r = 1:numel(roi)
    if makeplot==1
        figure('Position',[500 500 1200 300],'color','w');
    end
    for phase = 1:2 % make a sliding window looking for peaks within 5s batches
        start_ind = start_inds(phase);
        end_ind = end_inds(phase);
        
        % find peak in a specified window
        cdn = strcat([cond,'_plus_firstForages']);
        z = 1;
        for reg = 1:length(currModel) % iterate through regressors (reg = 1 for just doing PPI stats)
    %         if makeplot==1
    %             z = 1;
    %             figure('Position',[500 500 1200 300],'color','w');
    %         end
            all_means = [];
            all_peak_vals = [];
            all_sd = [];
            dfs = [];
            t = [];
            p = [];
%             all_p_wilcox = [];
%             all_zval_wilcox = [];
            avg_peak_sec = []; % keep track of when the peak is for each participant in seconds
            
            peak_indices = [];
            val_tzero = []; % keep track of what value the signal was at 0s for each participant, so we can use this for the t-test
            for s = 1:length(numsession) % iterate through participants
    
                numsession(s) = []; % drop this participant
                window  = mean(squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',cdn))(reg,:,numsession)),2); %avg across all timepoints for this participant/condition
    
                [m,i]  = max(abs(window(start_ind:end_ind)));
                i      = i + (start_ind-1);
                peak_indices = [peak_indices, i];
                val_tzero = [val_tzero, allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',cdn))(reg,22,s)]; %index 22 is the sample corresponding to zero
                peak(s,r) = allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',cdn))(reg,i,s);
    
                numsession = 1:numel(animals);
                clear window fw
    
            end
            avg_peak_sec = [avg_peak_sec, time(round(mean(peak_indices)))];
    
            [h,pv,ci,stats]= ttest(peak(:,r));
            all_peak_vals = [all_peak_vals, [peak(:,r)]];
            p = [p, pv];
            t = [t, stats.tstat];
            dfs = [dfs, stats.df];
            all_means = [all_means, mean(peak(:,r))];
            all_sd = [all_sd, stats.sd];
            
            % Also get non-parametric wilcoxon signed rank results
            [pv_wilcox,h_wilcox,stats_wilcox] = signrank(peak(:,r));
%             all_p_wilcox = [all_p_wilcox, pv_wilcox];
%             all_zval_wilcox = [all_zval_wilcox, stats_wilcox.zval];
            
            % Plot group aggregate with significance line
            if makeplot==1
                subplot(2, round(length(currModel)/2), z);
                beta = squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',cdn))(reg,:,:));
                stdshade(beta',0.15,plt_colour{1},linspace(-pre_win,post_win,nsamples)')
                if strcmp(currModel{reg}, 'seedTC')
                    ylim([-.2, .1])
                else
                    ylim([-0.1 0.1])
                end
                xlim([(-pre_win) post_win])
                set(gca,'fontsize',13)
                yL = get(gca,'YLim');
                line([0 0],yL,'Color',[0.8 0.8 0.8]);
                xL = get(gca,'XLim');
                line(xL,[0 0],'Color',[0.8 0.8 0.8]);

                if contains(regLabelMap(currModel{reg}),'seedTC')
                    title(strcat("Effect of ",replace(regLabelMap(currModel{reg}),'seedTC',roiLabelMap(seed_roi{r}))),'Interpreter','none','FontWeight','Normal');
                    subtitle(strcat("on ",roiLabelMap(roi{r})));
                else
                    title(strcat("Effect of ",regLabelMap(currModel{reg})),'Interpreter','none','FontWeight','Normal');
                    subtitle(strcat("on ",roiLabelMap(roi{r})));
                end

                xlabel(gca, 'Time (s)');
                ylabel(gca, 'Effect on BOLD signal');
                hold on
                %add stars for peaks added
%                 for pk = 1:length(numsession)
%                     clear plot
%                     plot(time(peak_indices(pk)),peak(pk,r),'b*');
%                     hold on
%                 end
                % Add star if this roi is significant (uncorrected)
                if pv_wilcox < 0.05
                    xline(time(round(mean(peak_indices))),'--p');
                    %plot(time(round(mean(peak_indices))),0.13,'b*');
                end
                z = z+1;
                ax = gca; 
                ax.FontSize = 16;
            end
            clear peak_indices val_tzero
            
            % False discovery rate
            [bonf_p, bonf_h] = bonf_holm(p);
%             [bonf_p_wilcox, bonf_h_wilcox] = bonf_holm(all_p_wilcox);
    
            % Log stats
            Condition = {cdn}; % repmat({cdn},1,length(roi))';
            Regressor = {currModel{reg}}; % repmat({currModel{reg}},1,length(roi))';
            ROI = {roi{r}};
            Seed = {seed_roi{r}};  % repmat({seed_roi{1}},1,length(roi))';
            Uncorrected_pvalue = p;
%             HolmBonf_corrected_pvalue = bonf_p';
            Significant = bonf_h;
            T_stats = t;
            df = dfs;
            Avg_peak_sec = avg_peak_sec;
            Peak_search_window = {strcat([num2str(round(time(start_ind),2)),'s-',num2str(round(time(end_ind),2)),'s'])};  % repmat({strcat([num2str(round(time(start_ind),2)),'s-',num2str(round(time(end_ind),2)),'s'])},1,length(roi))';
            Mean = all_means;
            SD = all_sd;
            p_wilcoxon = pv_wilcox;
%             p_wilcoxon_corrected = bonf_p_wilcox';
            zval_wilcoxon = stats_wilcox.zval;
            signedrank_wilcoxon = stats_wilcox.signedrank;
            Full_model = {strrep(strjoin(currModel),' ',' + ')};  % repmat({strrep(strjoin(currModel),' ',' + ')},1,length(roi))';
            Peak_vals = all_peak_vals';
    
            temp_stats = table(Condition,Regressor,ROI,Seed,Uncorrected_pvalue,Significant,T_stats,df,Avg_peak_sec,...
                Peak_search_window,Mean,SD,p_wilcoxon,zval_wilcoxon,signedrank_wilcoxon,Full_model,Peak_vals); % HolmBonf_corrected_pvalue, p_wilcoxon_uncorrected, p_wilcoxon_corrected, zval_wilcoxon
            allstats = [allstats;temp_stats];
            clear Condition Regressor ROI Seed Uncorrected_pvalue Significant temp_stats T_stats dfs Mean SD Peak_vals all_peak_vals signedrank_wilcoxon zval_wilcoxon p_wilcoxon% HolmBonf_corrected_pvalue
        end
    end
    if makeplot==1
        sgtitle({strcat("Condition: ",cdn,' | ROI: ',roi{r},' | Seed: ',seed_roi{r}),strcat('Early win: ',num2str(round(time(start_inds(1)),2)),'s-',num2str(round(time(end_inds(1)),2)),'s , Late win:',num2str(round(time(start_inds(2)),2)),'s-',num2str(round(time(end_inds(2)),2)),' | Test: Wilcoxon | Action coding: [-1,1] (normalized)')},'Interpreter','none');
        pos = get(gcf, 'Position');
        set(gcf, 'Position',pos+[500 0 500 200])
%         saveas(gcf,strcat('/Users/haileytrier/Desktop/PhD/foraging-under-threat/fMRI/Manuscript/ROI_and_PPI_analyses/Seed_',seed_roi{r},'_to_',roi{r},'_',regLabelMap(currModel{1}),'_',regLabelMap(currModel{2}),'_preWin',num2str(pre_win),'_NoLateWinTest'),'epsc');
    end
end
% savename=strcat('/Users/haileytrier/Desktop/PhD/foraging-under-threat/fMRI/Manuscript/ROI_and_PPI_analyses/PPI_',seed_roi,'x',currModel{2},'-to-',roi{r},'_checkVsForage_',cond,'.xlsx');
% writetable(allstats,savename{1},'WriteRowNames',false)

%% 4. Plot all time courses being tested for a particular regressor
% % First run sections 1 and 2.
% warning('off','all')
% 
% reg =  7; % index a regressor from currModel
% phase = 1; % index either early (1) or late (2) phase
% r = 1; % roi
% r_ind = '_1';
% disp(strcat("Displaying coefficients for reg = ",currModel{reg}," in phase ",num2str(phase)))
% 
% numsession = 1:numel(animals);
% post_win = 15;
% TR = 1.962;
% nsamples = round(((pre_win+post_win)./TR)*upsample);
% time = -pre_win:(pre_win+post_win)/nsamples:post_win;
% 
% if pre_win==2
%     start_inds = [22,72]; % if 2s prewindow: [22,72] for [0.0636s, 4.9769s]
%     end_inds = [72,123]; % if 2s prewindow: [72,123] for [4.9769s, 9.9884s]
% %     start_inds = [32,88]; % if 2s prewindow: [22,72] for [0.0636s, 4.9769s]
% %     end_inds = [88,124]; % if 2s prewindow: [72,144] for [4.9769s, 12.052]
% elseif pre_win==5
%     start_inds = [1,103]; % if 5s prewindow: [1,103] for -5-5s 
%     end_inds = [103,154]; % if 5s prewindow: [103,154] for 5-10s
% %     start_inds = [21, 73]; % [21, 73] for -3.0392-2.0588s
% %     end_inds = [73, 124]; % [73, 124] for 2.0588-7.0588s
% end
% 
% start_ind = start_inds(phase);
% end_ind = end_inds(phase);
% 
% % Set up plot
% clr = [0.8,0,0];
% plt_colour = {clr};
% z = 1;
% figure('Position',[500 500 1200 300],'color','w');
% 
% peak_indices = [];
% for s = 1:length(numsession) % iterate through participants
% 
%     numsession(s) = []; % drop this participant
%     window  = mean(squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},r_ind)).(sprintf('%s',strcat([cond,'_plus_firstForages'])))(reg,:,numsession)),2); %avg across all timepoints for this participant/condition
% 
%     % plot average of all other participants
%     subplot(5, 5, z); 
%     plotme = squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},r_ind)).(sprintf('%s',strcat([cond,'_plus_firstForages'])))(reg,:,numsession));
%     stdshade(plotme',0.15,plt_colour{1},linspace(-pre_win,post_win,nsamples)');
%     xlim([-pre_win post_win]);
%     hold on
% 
%     [m,i]  = max(abs(window(start_ind:end_ind)));
%     i      = i + (start_ind-1);
%     peak_indices = [peak_indices, i];
%     peak(s,r) = allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},r_ind)).(sprintf('%s',strcat([cond,'_plus_firstForages'])))(reg,i,s);
% 
%     % Plot point taken from this participant
%     plot(time(1:length(window)),allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},r_ind)).(sprintf('%s',strcat([cond,'_plus_firstForages'])))(reg,:,s),'Color','b');
%     plot(time(i),peak(s,r),'b*')
%     xline(time(i),'--p');
%     title(strcat("Participant ",num2str(s)," | Peak: ",num2str(peak(s,r))),'Interpreter','none');
% 
%     numsession = 1:numel(animals);
%     clear window fw
%     z = z+1;
% 
% end
% cdn = strcat([cond,'_plus_firstForages']);
% sgtitle(strcat("Condition: ",cdn,' | ROI: ',roi{r},' | Seed: ',seed_roi, " | Regressor: ",currModel{reg}),'Interpreter','none');
