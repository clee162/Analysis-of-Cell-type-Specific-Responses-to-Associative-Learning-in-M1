% for each ROI, find how many ITI lick bouts it responds to - then take average %. 
% Do the same for each reward 
% plot ITI response vs Reward response

clear
close all
clc

%% Data path 
save = 'n';
% 
celltype = 'PN';
mice_IDs ={'CL174', 'CL175', 'CL176', 'CL181', 'CL182', 'CL184'};

% celltype = 'PV';
% mice_IDs = {'CL147','CL172', 'CL195', 'CL196', 'CL198','CL158'};
% 
% celltype = 'VIP';
% mice_IDs = {'CL136', 'CL144', 'CL146', 'CL188'};

% celltype = 'SOM';
% mice_IDs = {'CL170', 'CL173', 'CL186', 'CL187', 'CL191', 'CL192', 'CL193'};

outer_folder = 'D:\2P data';
n = length(mice_IDs);

%% INPUTS
thresh_a = 1;       % active threshold
event_length = 5;   % frames of an active event 
t1_plot = -1;       % time in seconds before 1st lick for dF/F trace (for plot only)

%% initialize
n = n;
fr = 30;
window = 2.5;       % in seconds
days = {'Day 1', 'Day 7'};
num_days = length(days);
lick_bout_counter = zeros(n,length(days));
num_trials_all = zeros(n,length(days));
t_between_bouts = 3;

for d = 1:num_days
    disp(strcat('Analysing ', days{d}))
    for mouse = 1:n

        df = [];
        tone_start = [];
        wt_start = [];
        real_rois = [];
        num_reward_responses_daily = [];   %number of trials an ROI responds to on a given day. each row is an ROI
        num_ITI_responses_daily = [];

        %% load data
        curr_mouse = mice_IDs{mouse};
        curr_path = strcat(outer_folder,filesep,celltype,filesep,curr_mouse,filesep,days{d});
        disp(curr_path)
        matfiles =  dir(fullfile(curr_path, '*.mat'));
        nfiles = length(matfiles);
        for j = 1:length(matfiles)
            s = (fullfile(curr_path, matfiles(j).name));
            load(s)
        end

        %% exclude ROIs that don't fire at least once in the entire session
        if isempty(find(real_rois))
            disp([strcat(' mouse ', num2str(mouse), ' has no responsive rois')])
        else
            real_rois(real_rois == 0) = [];
            df = df(real_rois,:);
%             df(any(isnan(df),2),:) = [];
        end
        [num_rois, frames] = size(df);
        
        %% pre-processing 
        [df_z, wt_start, lick_start, tone_start] = daily_preprocessing(df, wt_start, lick_start, tone_start);
        num_trials = length(wt_start);
        num_trials_all(mouse,d) = length(wt_start);

        %% FIND REWARD RESPONSES
        % reward_responsive_ROI_indices is a cell array where each cell is
        % a trial
        [reward_responsive_ROI_indices, ~] = find_indices_of_reward_active_ROIs_v4(window, df_z, wt_start);
        reward_responsive_ROI_ind_mat = cell2mat(reward_responsive_ROI_indices); 
        reward_responsive_ROI_unique = unique(reward_responsive_ROI_ind_mat);
            
        %% FIND ITI LICKS
        ITI_responsive_ROI_indices = [];
        rois_abv_thresh = [];
        counter2 = 1;   %for lickbouts for each mouse. counter2 and lick_bout_counter are redundant. could make lick_bout_counter(mouse,d) = counter2
        for trial = 1:num_trials
            ITI_lick_ind = [];
            lick_bout_start_ind = [];
            perc_act_single_lickbout = [];
            ITI_responsive_ROI_ind_mat = [];
            ITI_active_roi = zeros(1,num_rois);
            
            if trial < num_trials    % if it's not the last trial
                ITI_lick_ind = find((lick_start > (wt_start(trial) + 20)) & (lick_start < (tone_start(trial+1)-window)));
            end
            if trial == num_trials    % if it's the last trial
                ITI_lick_ind = find((lick_start > (wt_start(trial) + 20)) & (lick_start < (frames/fr-window)));
            end
            
            % Find start of lick bouts
            curr_ITI_lick_times = lick_start(ITI_lick_ind);
            padded_curr_ITI_lick_times = vertcat(0, curr_ITI_lick_times, 0);
            diff_lick = diff(padded_curr_ITI_lick_times);
            lick_bout_start_ind = find(diff_lick > t_between_bouts);  %%
            
            if ~isempty(lick_bout_start_ind)
                 lick_bout_counter(mouse,d) = lick_bout_counter(mouse,d) + length(lick_bout_start_ind);

                 for j = 1:length(lick_bout_start_ind) 

                     ITI_lick_time = curr_ITI_lick_times(lick_bout_start_ind(j));
                     f1 = round(ITI_lick_time*fr);
                     f2 = f1+round(window*fr);
                     s_ITI = df_z(:,f1:f2);
                     
                     f1_plot_ITI(counter2) = f1 - fr;
                     f2_plot_ITI(counter2) = f2;
                     
                        for roi = 1:num_rois
                            frames_abv_thresh = zeros(length(real_rois),round(window*fr));  %initialize
                            frames_abv_thresh = (s_ITI(roi,:) > thresh_a);               % all frames above threshold
                            frames_abv_thresh = [false, frames_abv_thresh, false];
                            edges = diff(frames_abv_thresh); 
                            rising = find(edges==1);                                       % rising/falling edges
                            falling = find(edges==-1);  
                            spanwidth = falling - rising;                                  % how many frames between rising and falling edges
                            wideEnough = spanwidth >= event_length;                        % is the event long enough?
                            if ~isempty(find(wideEnough)) 
                                rois_abv_thresh(roi, j) = 1;

                            else 
                                rois_abv_thresh(roi,j) = 0;
                            end
                            if rois_abv_thresh(roi,j) > 0
                                ITI_active_roi(roi) = roi;
                            else 
                                ITI_active_roi(roi) = 0;
                            end
                        end
                     ITI_active_roi(ITI_active_roi == 0) = [];
                     ITI_responsive_ROI_indices{trial} = ITI_active_roi; 
                 end
            end
        end
        
        ITI_responsive_ROI_ind_mat = cell2mat(ITI_responsive_ROI_indices);
        ITI_responsive_ROI_ind_mat_unique = unique(ITI_responsive_ROI_ind_mat);
        
        for curr_roi = 1:num_rois
            num_reward_responses_daily(curr_roi) = length(find(reward_responsive_ROI_ind_mat == curr_roi));   %number of trials an ROI responds to on a given day. each row is an ROI
            num_ITI_responses_daily(curr_roi) = length(find(ITI_responsive_ROI_ind_mat == curr_roi));
        end
        
        num_reward_responses{mouse,d} = num_reward_responses_daily;
        num_ITI_responses{mouse,d} = num_ITI_responses_daily;
       
        percent_reward_responses(mouse,d) = mean((num_reward_responses_daily./num_trials_all(mouse,d)).*100);
        percent_ITI_responses(mouse,d) = mean((num_ITI_responses_daily./lick_bout_counter(mouse,d)).*100);  
        
    end
end

%% MEAN AND SEM OF % TRIALS/BOUTS 
SEM_ITI = (std(percent_ITI_responses,0,1))./sqrt(n);
SEM_reward = (std(percent_reward_responses,0,1))./sqrt(n);

overall_mean_ITI = mean(percent_ITI_responses,1);
overall_mean_reward = mean(percent_reward_responses,1);

%% make figures 
if strcmp(celltype, 'PN')
    c1 = [0.7 0.7 0.7];
    c2 = [0 0 0];
end
if strcmp(celltype, 'PV')
    c1 = [1 .5 .5];
    c2 = [.9 0 0];
end
if strcmp(celltype, 'VIP')
    c1 = [0.6 0.8 1];
    c2 = [0 0.5 1];
end
if strcmp(celltype, 'SOM')
    c1 = [120/255 220/255 120/255];
    c2 = [0/255 125/255 0/255];
end

x = [1 2];
y_star = 35;
figure(1)
for d = 1:num_days 
    subplot(1,num_days,d) 
    for mouse = 1:n 
        plot(x,[percent_ITI_responses(mouse,d) percent_reward_responses(mouse,d)],'-o',...
            'Color', c1,...
            'MarkerFaceColor', c1)
        hold on
    end
    e = errorbar(x, [overall_mean_ITI(d) overall_mean_reward(d)], [SEM_ITI(d) SEM_reward(d)], '-o', 'MarkerSize', 4,...
    'MarkerFaceColor', c2,... 
    'MarkerEdgeColor', c2,... 
    'LineWidth', 1.5);
    e.Color = c2;
    xlim([0 3]) 
    xticks([1 2])
    ylim([0 50])
    xticklabels({'ITI','reward'})
    title(days(d))
    [h,p(d)] = ttest(percent_ITI_responses(:,d), percent_reward_responses(:,d));
    NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
    text(NE(1), NE(2), ['p = ', num2str(p(d))], 'VerticalAlignment','top', 'HorizontalAlignment','right')
    if p(d) < 0.05
        text(max(xlim)/2.2, y_star,'*','FontSize',20)
    end 
    if d == 1
        ylabel('Avgerage ROI reliability (%)')
    end
end
suptitle(celltype)

%% save
file_name_e = strcat(celltype,'_ITI_vs_reward.eps'); % .eps OR .tif
file_name_t = strcat(celltype,'_ITI_vs_reward.tif'); % .eps OR .tif

if strcmp(save, 'Y') == 1
    save_path = 'C:\Users\clee162\Desktop\Figures\21-05-17 SOM figures\ITI SOM';
    s = strcat(save_path, filesep, file_name_e);
    print(gcf,'-depsc','-painters',s) % TO SAVE AS EPS
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
end

d1_I = 'day 1 ITI: ';
d7_I = 'day 7 ITI: ';

d1_R = 'day 1 reward: ';
d7_R = 'day 7 reward: ';
clc

%% display results
for d = 1:2
    
    if d == 1
        disp('Day 1')
        ITI_s = (strcat(d1_I, num2str(round(overall_mean_ITI(d),2)),{'+/-'}, num2str(round(SEM_ITI(d),2)),{'%, '}));
        R_s = strcat(d1_R, num2str(round(overall_mean_reward(d),2)),{'+/-'}, num2str(round(SEM_reward(d),2)),'%, ');
        p_s = strcat({'p = '}, num2str(round(p(d),2)));
        disp(strcat(ITI_s, R_s, p_s))
    end
    
    if d == 2
        disp('Day 7')
        ITI_s = (strcat(d7_I, num2str(round(overall_mean_ITI(d),2)),{'+/-'}, num2str(round(SEM_ITI(d),2)),{'%, '}));
        R_s = strcat(d7_R, num2str(round(overall_mean_reward(d),2)),{'+/-'}, num2str(round(SEM_reward(d),2)),'%');
        p_s = strcat({'p = '}, num2str(round(p(d),2)));
        disp(strcat(ITI_s, R_s, p_s))
    end
end
