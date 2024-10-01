%% 1. Select ROIs and model
clear all; fclose all; clc;
warning('on','all');
GLM = 3;
stat = 1;
textout = 0;
plot_session = 1;
plot_session_group = 0;
upsample = 20;
pre_win = 2;
data_root = '/Users/haileytrier/Downloads/Trier_et_al_2023_code/data/';
results_root = '/Users/haileytrier/Downloads/Trier_et_al_2023_code/results/';

timecourse_dir = strcat('/Users/haileytrier/Downloads/Trier_et_al_2023_code/data/timecourses/');
feat_name = 'L1stats_27Apr2022_noTempDerivs_motionRegressors_NoInterestContrasts';

subjects = {'sub103','sub105','sub109','sub111','sub117','sub119','sub121','sub123','sub125','sub127',...
    'sub129','sub131','sub204','sub208','sub212','sub214','sub216','sub218','sub222','sub224','sub228','sub230','sub232'};

roi = {'DRN_sub_PAG_edit_sub_SC_edit_ero','HB'};
seed_roi = {'none','none'};

% Load behavioural data
load(strcat(data_root,'all_subjects_behav_updatedEpochedData_upsample20_minusHalfTR_1Sep2024.mat'));

% Set model
currModel = {'switchToCheck','proximity','time','reward','noOfCones','constant','PPI1_noOfConesXswitchToCheck'};

% Map regressor names to nice labels for final figures
regnames = {'switchToCheck','timepressure','time','reward','seedTC','constant',...
    'PPI1_seedTCXswitchToCheck','PPI2_seedTCXtimepressure',...
    'PPI3_switchToCheckXtimepressure','PPI4_seedTCXswitchToCheckXtimepressure',...
    'switchToForage','PPI1_seedTCXswitchToForage','PPI2_seedTCXreward',...
    'PPI3_switchToForageXreward','PPI4_seedTCXswitchToForageXreward'...
    'PPI3_switchToForageXtimepressure','PPI4_seedTCXswitchToForageXtimepressure',...
    'PPI3_switchToCheckXreward','PPI4_seedTCXswitchToCheckXreward',...
    'PPI3_switchToCheckXproximity','PPI2_seedTCXproximity','PPI4_seedTCXswitchToCheckXproximity',...
    'proximity','noOfCones','PPI1_noOfConesXswitchToCheck'};
reglabels = {'Check','Time pressure','Time','Reward','seedTC','Constant',...
    'seedTC*Check','seedTC*Time pressure','Check*Time pressure','seedTC*Check*Time pressure',...
    'Forage','seedTC*Forage','seedTC*Reward','Forage*Reward','seedTC*Forage*Reward',...
    'Forage*Time pressure','seedTC*Forage*Time pressure','Check*Reward','seedTC*Check*Reward',...
    'Check*Proximity','seedTC*Proximity','seedTC*Check*Proximity','Proximity','Directions to Check','Check*Directions'};
regLabelMap = containers.Map(regnames, reglabels);
clear regnames, reglabels;

% Map ROI names to nice labels
roinames = {'striatal_sphere_rad4_dil','SN_Pauli',...
    'VTA_Nima_Pauli_sub_SN_Pauli','DRN_sub_PAG_edit_sub_SC_edit_ero',...
    'precentral_constant_rad2_7Feb','ACC_cluster_cope7','insular_cortex_cope7','HB'};
roilabels = {'Striatum','SN','VTA','DRN','Precentral gyrus','ACC','AI','Hb'};
roiLabelMap = containers.Map(roinames, roilabels);
clear roinames, roilabels;

%% 2. Run GLM
cond = 'monitoring_firstChecks_constant';

% Get data for this condition
firstChecksData = updatedBehav(strcmp(updatedBehav.action, 'firstChecks') & strcmp(updatedBehav.phase, 'monitoring') & updatedBehav.first_check==1, :);
firstForagesData = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'monitoring'), :);

for r = 1:length(roi)
    z=0;
    for a = 1:length(subjects)
        subject = subjects{a};
        z = z+1;
    
        subj_firstChecksData = firstChecksData(strcmp(firstChecksData.subjectNumber, subject), :);
        subj_firstForagesData = firstForagesData(strcmp(firstForagesData.subjectNumber, subject), :);
    
        % Run time course
        REG.switchToCheck = [repmat(1,[size(subj_firstChecksData,1),1]); repmat(-1,[size(subj_firstForagesData,1),1])]; %1=check, -1=forage
        REG.switchToForage = [repmat(-1,[size(subj_firstChecksData,1),1]); repmat(1,[size(subj_firstForagesData,1),1])]; %-1=check, 1=forage
        REG.constant = ones (length (REG.switchToCheck),1);
        REG.time = [subj_firstChecksData.('onset'); subj_firstForagesData.('onset')]; % checks, then forages
        REG.reward = [subj_firstChecksData.('reward'); subj_firstForagesData.('reward')]; % checks, then forages
        REG.timepressure = [subj_firstChecksData.('timePressure'); subj_firstForagesData.('timePressure')];
        REG.posUncertainty = [subj_firstChecksData.('posUncertainty'); subj_firstForagesData.('posUncertainty')];
        REG.proximity = [subj_firstChecksData.('proximity'); subj_firstForagesData.('proximity')];
        REG.noOfCones = [subj_firstChecksData.('noOfCones'); subj_firstForagesData.('noOfCones')];

        % load time-series data
        data_identifier = '_minusHalfTR_epoched';
        
        roi_TC_checks  = load([timecourse_dir,subject,'/',feat_name,'/',(['tc_',roi{r},'_',cond,'_upsample',num2str(upsample),data_identifier])]);
        roi_TC_forages  = load([timecourse_dir,subject,'/',feat_name,'/',(['tc_',roi{r},'_monitoring_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
        
        roi_TC = [roi_TC_checks.trial_data; roi_TC_forages.trial_data];

        clear roi_TC_checks roi_TC_forages
    
        % run GLM
        betaOut = NaN(size(currModel,2),size(roi_TC,2));
        for i = 1:size(roi_TC,2)
        
            % Make regressor design matrix (don't forget constant!)
            dmat = [];
            for c = 1:length(currModel)
                if ~contains(currModel{c},'PPI')
                    dmat = [dmat REG.(currModel{c})];
                end
            end
    
            % remove trials with no response
            dmat(isnan(dmat(:,2)),:)=[];
    
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1))); % normalize binary action var
            
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
            clear dmat contrasts REG.time REG.constant
        end
    
        allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',cond))(:,:,z) = betaOut;
        clear betaOut
    end
end
%% 3. Test, plot, and save stats
tableVarNames = {'Condition','Regressor','ROI','Seed','Uncorrected_pvalue','Significant',...
    'T_stats','df','Avg_peak_sec','Peak_search_window','Mean','SD','p_wilcoxon','signedrank_wilcoxon','Full_model',... % 'zval_wilcoxon'
    'Peak_vals'};
allstats = cell2table(cell(0,16), 'VariableNames', tableVarNames);
format shortG

clearvars window LOOT peak
numsession = 1:numel(subjects);
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
        cdn = cond;
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
    
                numsession = 1:numel(subjects);
                clear window fw
    
            end
            disp(strcat(roi{r},', phase=',num2str(phase),', reg=',currModel{reg}))
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
            
            % Plot group aggregate with significance line
            if makeplot==1
                subplot(2, round(length(currModel)/2), z);
                beta = squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',cdn))(reg,:,:));
                stdshade(beta',0.15,plt_colour{1},linspace(-pre_win,post_win,nsamples)')
                ylim([-.2, .2])
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
                % for pk = 1:length(numsession)
                %     clear plot
                %     plot(time(peak_indices(pk)),peak(pk,r),'b*');
                %     hold on
                % end
                % Add star if this roi is significant (uncorrected)
                if pv_wilcox < 0.05
                    xline(time(round(mean(peak_indices))),'--p');
                end
                z = z+1;
                ax = gca; 
                ax.FontSize = 16;
            end
            clear peak_indices val_tzero
            
            % False discovery rate
            [bonf_p, bonf_h] = bonf_holm(p);
    
            % Log stats
            Condition = {cdn};
            Regressor = {currModel{reg}}; 
            ROI = {roi{r}};
            Seed = {seed_roi{r}};
            Uncorrected_pvalue = p;
            Significant = bonf_h;
            T_stats = t;
            df = dfs;
            Avg_peak_sec = avg_peak_sec;
            Peak_search_window = {strcat([num2str(round(time(start_ind),2)),'s-',num2str(round(time(end_ind),2)),'s'])};  % repmat({strcat([num2str(round(time(start_ind),2)),'s-',num2str(round(time(end_ind),2)),'s'])},1,length(roi))';
            Mean = all_means;
            SD = all_sd;
            p_wilcoxon = pv_wilcox;
            %zval_wilcoxon = stats_wilcox.zval;
            signedrank_wilcoxon = stats_wilcox.signedrank;
            Full_model = {strrep(strjoin(currModel),' ',' + ')};  % repmat({strrep(strjoin(currModel),' ',' + ')},1,length(roi))';
            Peak_vals = all_peak_vals';
    
            temp_stats = table(Condition,Regressor,ROI,Seed,Uncorrected_pvalue,Significant,T_stats,df,Avg_peak_sec,...
                Peak_search_window,Mean,SD,p_wilcoxon,signedrank_wilcoxon,Full_model,Peak_vals); % HolmBonf_corrected_pvalue, p_wilcoxon_uncorrected, p_wilcoxon_corrected, zval_wilcoxon
            allstats = [allstats;temp_stats];
            clear Condition Regressor ROI Seed Uncorrected_pvalue Significant temp_stats T_stats dfs Mean SD Peak_vals all_peak_vals signedrank_wilcoxon zval_wilcoxon p_wilcoxon% HolmBonf_corrected_pvalue
        end
    end
    if makeplot==1
        sgtitle({strcat('Condition: Post-PD switch to check and switch to forage | ROI: ',roi{r},' | Seed: ',seed_roi{r}),strcat('Early win: ',num2str(round(time(start_inds(1)),2)),'s-',num2str(round(time(end_inds(1)),2)),'s , Late win:',num2str(round(time(start_inds(2)),2)),'s-',num2str(round(time(end_inds(2)),2)),' | Test: Wilcoxon | Action coding: [-1,1] (normalized)')},'Interpreter','none');
        pos = get(gcf, 'Position');
        set(gcf, 'Position',pos+[500 0 500 200])
    end
end