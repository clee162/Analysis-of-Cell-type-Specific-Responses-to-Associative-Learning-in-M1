% find the percent of trials each ROI responds to for tone and reward
% iterate through all days held in cell array "days" for a single mouse 
% INPUT: path = 'D:\2P data\PV\CL198' % character string; days =
% {'Day 1, 'Day 7'} % cell-array
% OUTPUT: percent_tone and percent_reward are matrices where the row index
% is the day specified by cell-array days and the column is the ROI index
% The value is the response % (or reliability) 
%-------------------------------------------------------------------------
% NOTES
% don't remove unresponsive ROIs from df matrix or else you mess up
% indexing. 
% replace unresponsive ROIs with NaN in the final output instead 


function [percent_tone,percent_reward] = updated_response_reliability(path,days)


fr = 30;
window = 2.5;       % in seconds
thresh_a = 1;       % active threshold
event_length = 5;   % frames of an active event 

for day = 1:length(days)
    curr_day = days(day);
    %% load data
    disp(strcat('Analysing   ',curr_day))

    df = [];
    tone_start = [];
    wt_start = [];
    real_rois = [];

    %% load data
    curr_path = strcat(path, filesep, curr_day);
    matfiles =  dir(fullfile(curr_path{1}, '*.mat'));
    nfiles = length(matfiles);
        for j = 1:length(matfiles)
            s = (fullfile(curr_path{1}, matfiles(j).name));
            load(s)
        end        
    unresponsive_ROIs = find(real_rois == 0);
    
    
    %% after day 1, remove any newly added ROIs
    % this is if you only want to analyze the ROIs that were there since
    % day 1
    [num_rois, frames] = size(df);
    if strcmp(curr_day, 'Day 1') == 1
        [day1_num_ROIs,~] = size(df);
    else
        df(day1_num_ROIs+1:num_rois,:) = [];
    end
    [num_rois, frames] = size(df);
    
    %% preprocessing and z-score
    [df_z, wt_start, lick_start, tone_start] = daily_preprocessing(df, wt_start, lick_start, tone_start);
    
    num_trials = length(tone_start);
        %% look for tone responsive ROIs
        for trial_i = 1:num_trials
            
            %tone frames
            f1_t = round(tone_start(trial_i)*fr);                         % tone onset frame
            f2_t = round(tone_start(trial_i)*fr)+round(window*fr)-1;      % tone+window frame
            df_z_t = df_z(:,f1_t:f2_t);
            
            % reward frames
            f1_r = round(wt_start(trial_i)*fr); 
            f2_r = f1_r + round(window*fr)-1;
            df_z_r = df_z(:,f1_r:f2_r);

            for roi = 1:num_rois

                % find tone responsive ROIs
                frames_abv_thresh_t = zeros(num_rois,round(window*fr)); %initialize
                frames_abv_thresh_t = (df_z_t(roi,:) > thresh_a);                  % all frames above threshold
                frames_abv_thresh_t = [false, frames_abv_thresh_t, false];
                edges_t = diff(frames_abv_thresh_t); 
                rising_t = find(edges_t==1);                                       % rising/falling edges
                falling_t = find(edges_t==-1);  
                spanwidth_t = falling_t - rising_t;                                % how many frames between rising and falling edges
                wideEnough_t = spanwidth_t >= event_length;                        % is the event long enough?
                if ~isempty(find(wideEnough_t)) 
                    tone_responsive_ROIs(roi, trial_i) = 1;
                else 
                    tone_responsive_ROIs(roi,trial_i) = 0;
                end


                % find reward responsive ROIs 
                frames_abv_thresh_r = zeros(num_rois,round(window*fr));            % initialize
                frames_abv_thresh_r = (df_z_r(roi,:) > thresh_a);                  % all frames above threshold
                frames_abv_thresh_r = [false, frames_abv_thresh_r, false];
                edges_r = diff(frames_abv_thresh_r); 
                rising_r = find(edges_r==1);                                       % rising/falling edges
                falling_r = find(edges_r==-1);  
                spanwidth_r = falling_r - rising_r;                                  % how many frames between rising and falling edges
                wideEnough_r = spanwidth_r >= event_length;                        % is the event long enough?
                if ~isempty(find(wideEnough_r)) 
                    reward_responsive_ROIs(roi, trial_i) = 1;
                else 
                    reward_responsive_ROIs(roi,trial_i) = 0;
                end

            end
        end


        %% To make the matrices the same size, add NaNs to previous day to fill difference in # of ROIs
%         if day > 1
%             if num_rois(day) > length(tone_trials) 
%                 new_rois = num_rois(day) - num_rois(day-1); 
%                 temp = NaN(day-1,new_rois);
%                 num_r_trials = horzcat(num_r_trials, temp);
%                 num_t_trials = horzcat(num_t_trials,temp);
%             end
%         end
        
   
        %% consolidate into output variables
        counter = 1;
        for roi = 1:day1_num_ROIs
            
            if counter > length(unresponsive_ROIs) || roi ~= unresponsive_ROIs(counter)
                num_t_trials(day,roi) = length(find(tone_responsive_ROIs(roi,:)));   %each element is a cell storing number of trials that ROI is active for
                num_r_trials(day,roi) = length(find(reward_responsive_ROIs(roi,:)));
                percent_tone(day,roi) = (num_t_trials(day,roi)./num_trials)*100;
                percent_reward(day,roi) = (num_r_trials(day,roi)./num_trials)*100;
            end
            
            % put NaN for unresponsive ROIs
            if counter <= length(unresponsive_ROIs)
                if roi == unresponsive_ROIs(counter)
                    num_t_trials(day,roi) = NaN;
                    num_r_trials(day,roi) = NaN;
                    percent_tone(day,roi) = NaN;
                    percent_reward(day,roi) = NaN;
                    counter = counter + 1;
                end
            end
            

        end
end



