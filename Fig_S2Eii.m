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

% Reward circuit: Replicate reward-related activity during switching
% to forage.
roi = {'thresh_zstat1_PrG','striatal_sphere_rad4_dil'};
model.pre = {'timepressure','time','reward','constant'};
model.post = {'proximity','time','reward','constant'};
myreg.constant = 4; % index of constant
myreg.parametric = 3; % index of threat
pre_win = 2;
GLM=2;
condition.pre = 'vigilance_firstForages_constant'; % loads first forage data from pre-disc. phase
condition.post = 'monitoring_firstForages_constant'; % loads first forage data from post-disc. phase
titles.constant = struct;
titles.parametric = struct;
titles.constant.figure = 'Ei. Switch to forage';
titles.parametric.figure = 'Eii. Reward during switch to forage';
titles.constant.filename = 'RewardCircuit-switchToForage-TC-replication_reproduced';
titles.parametric.filename = 'RewardCircuit-RewardDuringSwitchToForage-TC-replication_reproduced';
titles.statsfile = 'RewardCircuit-switchToForage-replication-stats_reproduced';
colors = 'greens';
ylims = [-0.4, 0.65];
searchWindow.start = [22,72]; % if 2s prewindow: [22,72] for [0.0636s, 4.9769s], [4.9769s, 9.9884s]
searchWindow.end = [72,123]; % if 2s prewindow: [72, 123] for [0.0636s, 4.9769s], [4.9769s, 9.9884s
twoTailPostTest = 0; % test post-disc. phase with a two-tailed test even if the pre-disc. phase is n.s.
subjects.post = {'sub105','sub109','sub111','sub119','sub121','sub125','sub127',...
    'sub131','sub216','sub224','sub230','sub232'}; % subjects with <2% effect required for reward in post-disc. checks & subjects with <40 checks omitted (N=12)

% Load concatenated behavioral data formatted for matlab
load(strcat(data_root,'all_subjects_behav_updatedEpochedData_upsample20_minusHalfTR.mat'));

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
    % Get data for this condition
    if strcmp(condition.(phase),'vigilance_firstChecks_constant')
        behavior = updatedBehav(strcmp(updatedBehav.action, 'allChecks') & strcmp(updatedBehav.phase, 'vigilance') & updatedBehav.first_check==1, :);
    
    elseif strcmp(condition.(phase),'monitoring_firstChecks_constant')
        behavior = updatedBehav(strcmp(updatedBehav.action, 'firstChecks') & strcmp(updatedBehav.phase, 'monitoring'), :);

    elseif strcmp(condition.(phase),'vigilance_firstForages_constant')
        behavior = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'vigilance'), :);
    
    elseif strcmp(condition.(phase),'monitoring_firstForages_constant')
        behavior = updatedBehav(strcmp(updatedBehav.action, 'firstForages') & strcmp(updatedBehav.phase, 'monitoring'),:);
    end
    
    for r = 1:length(roi)
        z=0;
        for a = 1:length(subjects.(phase))
            subj = subjects.(phase){a};
            z = z+1;
        
            subj_behavior = behavior(strcmp(behavior.subjectNumber, subj), :);

            % Run time course
            REG.constant = ones (length (subj_behavior.('onset')),1);
            REG.time = [subj_behavior.('onset')];
            REG.reward = [subj_behavior.('reward')];
            REG.timepressure = [subj_behavior.('timePressure')];
            REG.posUncertainty = [subj_behavior.('posUncertainty')];
            REG.proximity = [subj_behavior.('proximity')];
    
            % load time-series data
            if pre_win==5 data_identifier='_minusHalfTR_5sPreWindow_epoched';, else data_identifier='_minusHalfTR_epoched';, end            
            roi_TC  = load([timecourse_dir,subj,'/',feat_name,'/',(['tc_',roi{r},'_',condition.(phase),'_upsample',num2str(upsample),data_identifier])]);
            roi_TC = [roi_TC.trial_data];

            % run GLM
            betaOut = NaN(size(model.(phase),2),size(roi_TC,2));
            for i = 1:size(roi_TC,2)

                % Make regressor design matrix (don't forget constant!)
                dmat = [];
                for c = 1:length(model.(phase))
                    if ~contains(model.(phase){c},'PPI')
                        dmat = [dmat REG.(model.(phase){c})];
                    end
                end
        
                % remove trials with no response
                [nanrows,nancols]=find(isnan(dmat));
                dmat(nanrows,:) = [];
                clear nanrows nancols
        
                % normalise data
                dmat(:,1:myreg.constant-1) = zscore(dmat(:,1:myreg.constant-1)); % exclude constant from normalization
                dmat(:,myreg.constant+1:end) = zscore(dmat(:,myreg.constant+1:end));
                
                % Constrain to later post-disc. checks - msec_since_disc
%                     if strcmp(phase,'post')
%                         dmat = dmat(ind_top90,:);
%                     end

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));
        
                % beta X time output
                betaOut(:,i) = ols(roi_TC(:,i),dmat,contrasts);
                clear dmat contrasts
            end
        
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r},'_',num2str(r))).(sprintf('%s',condition.(phase)))(:,:,z) = betaOut;
            clear betaOut REG
        end
    end
    clear subj c z a r i 
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

%% Plot nested bars
% close all
fontsize = 25;
peaks = struct;
peaks.pre = struct;
peaks.post = struct;
for r = 1:numel(roi)
    peaks.pre.(roi{r}) = struct;
    peaks.post.(roi{r}) = struct;

    % Extract stats from statistics table
    for regType = {'constant', 'parametric'}
        regType = regType{1};
        rowindex = find(contains(allstats.Condition, condition.pre) & strcmp(allstats.Regressor, model.pre{myreg.(regType)}) & strcmp(allstats.ROI, roi{r}));
        peaks.pre.(roi{r}).(regType) = struct;
        peaks.pre.(roi{r}).(regType).values = table2array(allstats(rowindex,size(allstats,2)));
        peaks.pre.(roi{r}).(regType).p = table2array(allstats(rowindex,'p_wilcox')); % p value from two-sided wilcoxon test

        rowindex = find(contains(allstats.Condition, condition.post) & strcmp(allstats.Regressor, model.post{myreg.(regType)}) & strcmp(allstats.ROI, roi{r}));
        peaks.post.(roi{r}).(regType) = struct;
        peaks.post.(roi{r}).(regType).values = table2array(allstats(rowindex,size(allstats,2)));

        % Get p values for post-disc. stats according to sign of pre-disc. stats
        for peak = 1:numel(peaks.post.(roi{r}).(regType).values)
            if peaks.pre.(roi{r}).(regType).p{1} < 0.05
                peaks.post.(roi{r}).(regType).p = table2array(allstats(rowindex,'p_wilcox'));
            else
                peaks.post.(roi{r}).(regType).p = 'NA'; % don't report post-disc. test if pre-disc. test not significant
            end
        end
    end
end

% Plot
% Though we have calculated stats for the constant and parameteric terms, we're only
% plotting the parametric term here. For this panel we used the ROI mask corresponding 
% to reward-related activity in the precentral gyrus, 'thresh_zstat1_PrG', in 
% order to replicate reward-related results. For panel 7Bi we used the
% general precentral gyrus mask. See Fig_7Bi.m.
for targetReg = {'parametric'}
    targetReg = targetReg{1};
    figure

    % Get pre value stats
    prevals.mean = [];
    prevals.sem = [];
    for r = 1:numel(roi)
        prevals.mean = [prevals.mean; mean(peaks.pre.(roi{r}).(targetReg).values,'omitnan')];
        prevals.sem = [prevals.sem; std(peaks.pre.(roi{r}).(targetReg).values,'omitnan')/sqrt(length(peaks.pre.(roi{r}).(targetReg).values))];
    end

    % Get post value stats
    postvals.mean = [];
    postvals.sem = [];
    for r = 1:numel(roi)
        postvals.mean = [postvals.mean; mean(peaks.post.(roi{r}).(targetReg).values,'omitnan')];
        postvals.sem = [postvals.sem; std(peaks.post.(roi{r}).(targetReg).values,'omitnan')/sqrt(length(peaks.post.(roi{r}).(targetReg).values))];
    end

    vals.mean = [prevals.mean postvals.mean];
    vals.sem = [prevals.sem postvals.sem];
    clear prevals postvals

    h=bar(vals.mean);
    hold on
    errorbar(h(1).XEndPoints,vals.mean(:,1),vals.sem(:,1),'LineStyle','none','Color','k','LineWidth',2);
    errorbar(h(2).XEndPoints,vals.mean(:,2),vals.sem(:,2),'LineStyle','none','Color','k','LineWidth',2);
    if strcmp(colors,'oranges') h(1).FaceColor='#fac205';, else h(1).FaceColor='#acbf69';, end
    if strcmp(colors,'oranges') h(2).FaceColor='#c04e01';, else h(2).FaceColor='#028f1e';, end

    % Significance stars (wilcoxon tests)
    for r = 1:numel(roi)
        if vals.mean(r,1)<0 starheight=.1;, else starheight=vals.mean(r,1)+vals.sem(r,1)+0.03;,end
        text(r-0.15, starheight, calcSigStars(peaks.pre.(roi{r}).(targetReg).p{1}),'FontSize',20,'HorizontalAlignment','center');
        if peaks.pre.(roi{r}).(targetReg).p{1}<0.05
            if vals.mean(r,2)<0 starheight=.1;, else starheight=vals.mean(r,2)+vals.sem(r,2)+.03;,end
            text(r+0.15, starheight, calcSigStars(peaks.post.(roi{r}).(targetReg).p{1}),'FontSize',20,'HorizontalAlignment','center');
        end
        if strcmp(colors,'oranges') mfc={'#bd8e0d','#ff8000'};, else mfc={'#d0ff00','#25e000'};, end
        s1=scatter(r-0.145, peaks.pre.(roi{r}).(targetReg).values, 'MarkerFaceColor',mfc{1},'MarkerEdgeColor','#41423d','MarkerFaceAlpha',0.6,'jitter','on', 'jitterAmount', 0.05);
        s2=scatter(r+0.145, peaks.post.(roi{r}).(targetReg).values, 'MarkerFaceColor',mfc{2},'MarkerEdgeColor','#41423d','MarkerFaceAlpha',0.6,'jitter','on', 'jitterAmount', 0.05);
    end

    % Set other plot attributes
    title(titles.(targetReg).figure, 'FontSize', fontsize);
    xticklabels(cellfun(@(x) roilabel.(x), roi, 'UniformOutput', false))
    xlabel("ROI", 'FontSize', fontsize)
    ylabel("Peak values", 'FontSize', fontsize)
    ax = gca; 
    ax.XAxis.FontSize = fontsize;
    ax.YAxis.FontSize = fontsize;
%         yl = ylim;
%         ylim([min(min(vals.mean))-0.07, max(max(vals.mean))+0.07]); % make sure plot will show highest asterisk
    ylim(ylims);
end