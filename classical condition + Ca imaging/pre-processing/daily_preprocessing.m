function [df_z, wt_start, lick_start, tone_start] = daily_preprocessing(df, wt_start, lick_start, tone_start);

[num_rois, frames] = size(df);
fr = 30;

%% remove tone_starts that are past the imaging
tone_start = tone_start(tone_start < (frames/fr - 5*fr));
wt_start = wt_start(1:length(tone_start));

%% exclude trials with no lick response
[wt_start, lick_start, tone_start, num_trials] = find_trials_with_lick_response(wt_start, lick_start, tone_start);            

%% Get baseline and values for z-score 
[mu, sigma, baseline] = baseline_calculation(df, tone_start);

%% z-score 
df_z = zeros(num_rois, frames);
for roi = 1:num_rois
    for f = 1:frames
        df_z(roi,f) = (df(roi,f)-mu(roi))./sigma(roi);
    end
end