%% 1. Select ROIs and model

clear all; close all; clc;
warning('on','all');
upsample = 20;
data_root = '/Users/haileytrier/Downloads/Trier_et_al_2023_code/data/';
results_root = '/Users/haileytrier/Downloads/Trier_et_al_2023_code/results/';

timecourse_dir = strcat('/Users/haileytrier/Downloads/Trier_et_al_2023_code/data/timecourses/');
feat_name = 'L1stats_27Apr2022_noTempDerivs_motionRegressors_NoInterestContrasts';

subjects.pre = {'sub103','sub105','sub109','sub111','sub117','sub119','sub121','sub123','sub125','sub127',...
    'sub129','sub131','sub204','sub208','sub212','sub214','sub216','sub218','sub222','sub224','sub228','sub230','sub232'};

% Fig. 3C: Switch to check: Effect of ACC on Hb
roi = {'HB'};
seed_roi = {'ACC_cluster_cope7'};
titles.figuretitle = 'Switch to check: Effect of ACC on Hb';
titles.figurefile = 'Fig7C_ACC_to_Hb_during_checks-TC-replication';
titles.statsfile = 'Fig7C_ACC_to_Hb_during_checks-TC-replication-stats_reproduced';
model.pre = {'timepressure','time','reward','seedTC','constant','PPI2_seedTCXtimepressure'};
model.post = {'proximity','time','reward','seedTC','constant','PPI2_seedTCXproximity'};
myreg.constant = 5; % index of constant
myreg.parametric = 4; % index of parametric regressor of interest
pre_win = 5;
condition.pre = 'vigilance'; % loads first check and first forage data from pre-disc. phase
condition.post = 'monitoring'; % loads first check and first forage data from post-disc. phase
searchWindow.start = [1]; % if 5s prewindow: [1,103] for -5-5s 
searchWindow.end = [103]; % if 5s prewindow: [103,154] for 5-10s
twoTailPostTest = 0; % test post-disc. phase with a two-tailed test even if the pre-disc. phase is n.s.
subjects.post = {'sub105','sub109','sub111','sub119','sub121','sub123','sub125','sub127',...
    'sub131','sub216','sub224','sub230','sub232'}; % subjects with <2% effect required for constant in post-disc. checks & subjects with <40 checks omitted (N=13)
GLM=4;
condition.action = 'checks'; % which action transition to analyze
colors = 'oranges';

% Load concatenated behavioral data formatted for matlab
load(strcat(data_root,'all_subjects_behav_updatedEpochedData_upsample20_minusHalfTR_5sPreWindow.mat'));

% Map regressor names to nice labels for figures
[regLabel.switchToCheck, regLabel.timepressure, regLabel.time] = deal('Check','Threat level', 'Time');
[regLabel.reward, regLabel.seedTC, regLabel.constant] = deal('Reward','seedTC','Constant');
[regLabel.PPI1_seedTCXswitchToCheck, regLabel.PPI2_seedTCXtimepressure] = deal('seedTC*Check','seedTC*Threat level');
[regLabel.PPI3_switchToCheckXtimepressure, regLabel.PPI4_seedTCXswitchToCheckXtimepressure] = deal('Check*Threat level','seedTC*Check*Threat level');
[regLabel.switchToForage, regLabel.PPI1_seedTCXswitchToForage, regLabel.PPI2_seedTCXreward] = deal('Forage','seedTC*Forage','seedTC*Reward');
[regLabel.PPI3_switchToForageXreward, regLabel.PPI4_seedTCXswitchToForageXreward] = deal('Forage*Reward','seedTC*Forage*Reward');
[regLabel.PPI3_switchToForageXtimepressure, regLabel.PPI4_seedTCXswitchToForageXtimepressure] = deal('Forage*Threat level','seedTC*Forage*Threat level');
[regLabel.PPI3_switchToCheckXreward, regLabel.PPI4_seedTCXswitchToCheckXreward] = deal('Check*Reward','seedTC*Check*Reward');
[regLabel.PPI3_switchToCheckXproximity, regLabel.proximity, regLabel.posUncertainty] = deal('Check*Threat level','Threat level','Threat level');
[regLabel.PPI4_seedTCXswitchToCheckXproximity, regLabel.PPI4_seedTCXswitchToCheckXposUncertainty] = deal('seedTC*Check*Threat level','seedTC*Check*Threat level');
[regLabel.PPI2_seedTCXposUncertainty, regLabel.PPI2_seedTCXproximity] = deal('seedTC*Threat level','seedTC*Threat level');

% Map ROI names to nice labels
[roilabel.striatal_sphere_rad4_dil, roilabel.SN_Pauli, roilabel.VTA_Nima_Pauli_sub_SN_Pauli] = deal('Striatum','SN','VTA');
[roilabel.DRN_sub_PAG_edit_sub_SC_edit_ero, roilabel.precentral_constant_rad2_7Feb, roilabel.ACC_cluster_cope7] = deal('DRN','Precentral gyrus','ACC');
[roilabel.insular_cortex_cope7, roilabel.HB, roilabel.thresh_zstat1_PrG] = deal('AI','Hb','Precentral gyrus');

%% GLM
for phase = {'pre','post'} % iterate through pre- and post-discovery phases
    phase = phase{1};

    % Get behavioral data
    if strcmp(condition.action, 'checks')
        if strcmp(phase,'pre')
            all_behavior = updatedBehav(strcmp(updatedBehav.action, 'allChecks') & strcmp(updatedBehav.phase, 'vigilance') & updatedBehav.first_check==1, :);
        elseif strcmp(phase,'post')
            all_behavior = updatedBehav(strcmp(updatedBehav.action, 'firstChecks') & strcmp(updatedBehav.phase, 'monitoring'), :);
        end

    elseif strcmp(condition.action, 'forages')
        all_behavior = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, condition.(phase)), :);
    end

    for r = 1:length(roi)  % r=1;z=0;a=1;
        z=0;
        for a = 1:length(subjects.(phase))
            subj = subjects.(phase){a};
            z = z+1;
        
            subj_behavior = all_behavior(strcmp(all_behavior.subjectNumber, subj), :);

            % Run time course
            REG.constant = ones (length (subj_behavior.('onset')),1);
            REG.time = [subj_behavior.('onset')];
            REG.reward = [subj_behavior.('reward')];
            REG.timepressure = [subj_behavior.('timePressure')];
            REG.posUncertainty = [subj_behavior.('posUncertainty')];
            REG.proximity = [subj_behavior.('proximity')];
    
            % load time-series data
            if pre_win==5 data_identifier='_minusHalfTR_5sPreWindow_epoched';, else data_identifier='_minusHalfTR_epoched';,end
        
            if strcmp(condition.action, 'checks')
                seed_TC = load([timecourse_dir,subj,'/',feat_name,'/',(['tc_',seed_roi{r},'_',condition.(phase),'_firstChecks_constant_upsample',num2str(upsample),data_identifier])]);
                roi_TC  = load([timecourse_dir,subj,'/',feat_name,'/',(['tc_',roi{r},'_',condition.(phase),'_firstChecks_constant_upsample',num2str(upsample),data_identifier])]);
            elseif strcmp(condition.action, 'forages')
                seed_TC = load([timecourse_dir,subj,'/',feat_name,'/',(['tc_',seed_roi{r},'_',condition.(phase),'_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
                roi_TC  = load([timecourse_dir,subj,'/',feat_name,'/',(['tc_',roi{r},'_',condition.(phase),'_firstForages_constant_upsample',num2str(upsample),data_identifier])]);
            end
            roi_TC = [roi_TC.trial_data];
            seed_TC = [seed_TC.trial_data];
        
            % run GLM
            betaOut = NaN(size(model.(phase),2),size(roi_TC,2));
            for i = 1:size(roi_TC,2)
        
                % load ROI time-series data
                REG.seedTC  = seed_TC(:,i);

                % Make regressor design matrix (don't forget constant!)
                dmat = [];
                for c = 1:length(model.(phase))
                    if ~contains(model.(phase){c},'PPI')
                        dmat = [dmat REG.(model.(phase){c})];
                    end
                end
        
                % remove trials with no response
                dmat(isnan(dmat(:,2)),:)=[];
        
                % normalise data
                dmat(:,1:myreg.constant-1) = zscore(dmat(:,1:myreg.constant-1)); % exclude constant from normalization
                dmat(:,myreg.constant+1:end) = zscore(dmat(:,myreg.constant+1:end));
                
                % Add PPI regressors
                for ppi = 1:length(model.(phase))
                    if contains(model.(phase){ppi},'PPI')
                        % Determine vars that go into this PPI interaction
                        ppivars = split(model.(phase){ppi},'_');
                        ppivars = split(ppivars{2},'X');
                        if length(ppivars)==2
                            REG.(model.(phase){ppi}) = zscore(REG.(ppivars{1}) .* REG.(ppivars{2}));
                        elseif length(ppivars)==3
                            REG.(model.(phase){ppi}) = zscore(REG.(ppivars{1}) .* REG.(ppivars{2}) .* REG.(ppivars{3}));
                        end
                        dmat = [dmat, REG.(model.(phase){ppi})];
                    end
                end
                
                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));
        
                % beta X time output
                betaOut(:,i) = ols(roi_TC(:,i),dmat,contrasts);
                clear dmat contrasts
            end
        
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',strcat([condition.(phase)])))(:,:,z) = betaOut;
            clear betaOut REG
        end
    end
    clear subj c z a r all_behavior i roi_TC seed_TC subj_behavior data_identifier
end

%% Test and save stats
allstats = table; 
format shortG

for phase = {'pre','post'} % iterate through pre- and post-discovery phases
    phase = phase{1};
    clearvars window LOOT peak

    numsession = 1:numel(subjects.(phase));
    post_win = 15;
    TR = 1.962;
    nsamples = round(((pre_win+post_win)./TR)*upsample);
    time = -pre_win:(pre_win+post_win)/nsamples:post_win;
    
    for r = 1:numel(roi)
        if strcmp(phase,'pre') prephase.(roi{r}) = struct;,end

        % find peak in a specified window
        z = 1;
        for reg = [myreg.constant, myreg.parametric] % Target regressor(s)
            if strcmp(phase,'pre') prephase.(roi{r}).(strcat('reg',num2str(reg))) = struct;,end

            for win = 1 % make a sliding window looking for peaks within 5s batches
                if strcmp(phase,'pre') prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))) = struct;,end

                % Set indices for search window
                if strcmp(phase,'pre')
                    start_ind = searchWindow.start(win);
                    end_ind = searchWindow.end(win);
                elseif strcmp(phase,'post')
                    if twoTailPostTest
                        start_ind = searchWindow.start(2); % for fig. 6Aiii
                        end_ind = searchWindow.end(2); % for fig. 6Aiii
                    else
                        % Get the range within which 95% of pre-disc. phase peak indices lie
                        dist = abs(prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).peak_indices - mean(prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).peak_indices));
                        [sortDist, sortIndex] = sort(dist);
                        index_95perc = sortIndex(1:floor(0.95 * numel(prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).peak_indices)));
                        x_95percent = prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).peak_indices(index_95perc);
                        start_ind = min(x_95percent);
                        end_ind = max(x_95percent);
                    end
                end
                
                all_means = [];
                all_peak_vals = [];
                all_sd = [];
                dfs = [];
                t = [];
                p = [];
                avg_peak_sec = []; % keep track of when the peak is for each participant in seconds
    
                peak_indices = [];
                for s = 1:length(numsession) % iterate through participants
        
                    numsession(s) = []; % drop this participant
                    window  = mean(squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',condition.(phase)))(reg,:,numsession)),2);
                    
                    % Search for peak
                    [m,i]  = max(abs(window(start_ind:end_ind)));
                    
                    i      = i + (start_ind-1);
                    peak_indices = [peak_indices, i];
                    peak(s,r) = allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',condition.(phase)))(reg,i,s);
        
                    numsession = 1:numel(subjects.(phase));
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
                
                % Get non-parametric wilcoxon signed rank results
                if kstest(peak(:,r))==1 normalDist='False';, else normalDist='True';, end
                if strcmp(phase,'pre')
                    [pv_wilcox,h_wilcox,stats_wilcox] = signrank(peak(:,r),0);
                    prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).pval = pv_wilcox;
                    prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).tail = mean(peak(:,r));
                    prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).peak_indices = peak_indices;
                    tail='two_sided';
                elseif strcmp(phase,'post')
                    if prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).pval<0.05
                        if prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).tail>0
                            [pv_wilcox,h_wilcox,stats_wilcox] = signrank(peak(:,r),0,'tail','right');
                            tail='right';
                        elseif prephase.(roi{r}).(strcat('reg',num2str(reg))).(strcat('win_',num2str(win))).tail<0
                            [pv_wilcox,h_wilcox,stats_wilcox] = signrank(peak(:,r),0,'tail','left');
                            tail='left';
                        end
                    else
                        if twoTailPostTest
                            [pv_wilcox,h_wilcox,stats_wilcox] = signrank(peak(:,r),0);
                            tail='two-sided';
                        else
                            pv_wilcox='NA';
                            stats_wilcox.signedrank='NA';
                            stats_wilcox.zval='NA';
                        end
                    end
                end
                                    
                % False discovery rate
                [bonf_p, bonf_h] = bonf_holm(p);
        
                % Log stats
                temp_stats = table; 
                temp_stats.Condition = {condition.(phase)};
                temp_stats.Regressor = {model.(phase){reg}};
                temp_stats.ROI = {roi{r}};
                if exist('seed_roi','var') temp_stats.Seed={seed_roi{r}};, else tempm_stats.Seed={'NA'};, end;
                temp_stats.normalDist = {normalDist};
                temp_stats.N = {numel(all_peak_vals)};
                temp_stats.Uncorrected_pvalue = {p};
                temp_stats.Significant = {bonf_h};
                temp_stats.T_stats = {t};
                temp_stats.df = {dfs};
                temp_stats.Avg_peak_sec = {avg_peak_sec};
                temp_stats.Peak_search_window = {strcat([num2str(round(time(start_ind),2)),'s-',num2str(round(time(end_ind),2)),'s'])};  % repmat({strcat([num2str(round(time(start_ind),2)),'s-',num2str(round(time(end_ind),2)),'s'])},1,length(roi))';
                temp_stats.Mean = {all_means};
                temp_stats.SD = {all_sd};
                temp_stats.p_wilcox = {pv_wilcox};
                if length(fieldnames(stats_wilcox))==2 temp_stats.zval_wilcox={stats_wilcox.zval};, else temp_stats.zval_wilcox={'NA'};,end; % zval not returned when n<16
                temp_stats.signedrank_wilcox = {stats_wilcox.signedrank};
                temp_stats.tail_wilcox={tail};
                temp_stats.Full_model = {strrep(strjoin(model.(phase)),' ',' + ')};
                if length(all_peak_vals)<23
                    all_peak_vals(end+1:23) = missing;
                end
                temp_stats.Peak_vals = all_peak_vals';
                allstats = [allstats; temp_stats];
            end
        end
    end
end
savename=strcat(results_root,titles.statsfile,'.xlsx');
% writetable(allstats,savename,'WriteRowNames',false);

%% Plot bars
% close all
fontsize = 25;
stars_fontsize = 30;
prePeaks = table2array(allstats(strcmp(allstats.Condition, 'vigilance') & strcmp(allstats.Regressor,model.pre{myreg.parametric}),'Peak_vals'),size(allstats,2));
pval_pre = table2array(allstats(strcmp(allstats.Condition, 'vigilance') & strcmp(allstats.Regressor,model.pre{myreg.parametric}),'p_wilcox'));
postPeaks = table2array(allstats(strcmp(allstats.Condition, 'monitoring') & strcmp(allstats.Regressor,model.post{myreg.parametric}),size(allstats,2)));
pval_post = table2array(allstats(strcmp(allstats.Condition, 'monitoring') & strcmp(allstats.Regressor,model.post{myreg.parametric}),'p_wilcox'));

figure
h=bar([mean(prePeaks,'omitnan'); mean(postPeaks,'omitnan')],'FaceColor','flat');
if strcmp(colors,'oranges') h.CData(2,:) = [0.8500 0.3250 0.0980];, else h.CData(2,:)=[0.0078125, 0.55859, 0.11719];, end
if strcmp(colors,'oranges') h.CData(1,:) = [0.9290 0.6940 0.1250];, else h.CData(1,:)=[0.67188, 0.74609, 0.41016];, end
hold on;

errorbar(h(1).XEndPoints,[mean(prePeaks,'omitnan'); mean(postPeaks,'omitnan')], ...
    [std(prePeaks,'omitnan')/sqrt(length(prePeaks)); std(postPeaks,'omitnan')/sqrt(length(postPeaks))], ...
    'LineStyle','none','Color','k','LineWidth',2)

% Scatter point with jitter
if strcmp(colors,'oranges') preColor="#bd8e0d";, else preColor='#d0ff00';, end
if strcmp(colors,'oranges') postColor="#ff8000";, else postColor='#25e000';, end
scatter(repmat(1, length(prePeaks), 1),prePeaks,60,'MarkerFaceColor',preColor,'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.05,'MarkerFaceAlpha',0.6,'jitter','on', 'jitterAmount', 0.1) % scatter pre peaks
scatter(repmat(2, length(postPeaks), 1),postPeaks,60,'MarkerFaceColor',postColor,'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.05,'MarkerFaceAlpha',0.6,'jitter','on', 'jitterAmount', 0.1) % scatter post peaks

% Significance stars (wilcoxon tests)
text(1, max(prePeaks)+0.05, calcSigStars(pval_pre{1}),'FontSize',stars_fontsize,'HorizontalAlignment','center');
text(2, max(postPeaks)+0.05, calcSigStars(pval_post{1}),'FontSize',stars_fontsize,'HorizontalAlignment','center');

% Set other plot attributes
title(titles.figuretitle, 'FontSize', fontsize);
xticklabels(["Pre" "Post"])
xlabel("Discovery phase", 'FontSize', fontsize)
ylabel("Peak values", 'FontSize', fontsize)
ax = gca; 
ax.XAxis.FontSize = fontsize;
ax.YAxis.FontSize = fontsize;
yl = ylim;
ylim([yl(1), max([max(prePeaks)+0.1, max(postPeaks)+0.2])]); % make sure plot will show highest asterisk