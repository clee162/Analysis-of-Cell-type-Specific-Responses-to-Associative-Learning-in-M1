% (1) concatenate all baseline periods (2.5 s preceeding the tone for each trial) 
% (2) calculate mu and sigma to use for z-score 
% (3) average the concatenated baselines to obtain a baseline matrix of
% dimensions [Number of ROIs x window*fr]
% Based off Kato et al., Neuron, 2015


function [mu, sigma, baseline] = baseline_calculation(df, tone_start)

 %% calculate baseline mu and sigma for each ROI by horizontally concatenating all baselines
 
 %initialize
[num_rois,num_frames] = size(df);
mu = [];    % 1 x number of ROIs vector with each ROI's mean df/f during the baseline period
sigma = []; % 1 x number of ROIs vector with each ROI's standard deviation of df/f during baseline
baseline_cat_traces = [];
window = 2.5;   % window s before tone up to tone onset = baseline period 
fr = 30; 
session_length_s = num_frames/fr;
tone_start = tone_start(tone_start < (session_length_s - window)); %don't use tone_starts longer than imaging session
num_trials = length(tone_start);

 for trial = 1:num_trials
     f1_base = (round((tone_start(trial) - window) * fr));       % baseline start is 'window' s before tone
     f2_base = f1_base + round(window*fr) - 1;                   % baseline ends 1 frame before tone
     single_trial_baseline = df(:,f1_base:f2_base);
     
     if trial == 1
         baseline_cat_traces = single_trial_baseline;
     else
         baseline_cat_traces = horzcat(baseline_cat_traces, single_trial_baseline);
     end
     
 end
 sigma = std(baseline_cat_traces,0,2);                           % std of all baselines
 mu(:,1) = mean(baseline_cat_traces,2);                          % find a baseline mu value for each ROI
 
 %% calculate a baseline trace for each ROI by averaging it's baseline (frame by frame) for each trial
 baseline = [];
 for roi = 1:num_rois
     single_roi_all_baselines = [];
     for trial = 1:num_trials
         f1_base = (round((tone_start(trial) - window) * fr));
         f2_base = f1_base + round(window*fr) - 1;
         single_roi_baseline = df(roi,f1_base:f2_base);
         if trial == 1
             single_roi_all_baselines = single_roi_baseline;
         else
             single_roi_all_baselines = vertcat(single_roi_all_baselines, single_roi_baseline);
         end
     end
     baseline(roi,:) = mean(single_roi_all_baselines,1);
 end
 