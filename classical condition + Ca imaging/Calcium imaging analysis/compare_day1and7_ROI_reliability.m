% INPUT:  celltype = string 'VIP', 'PV', 'SOM' or 'PN';
% mice_IDs = cell array with different mice IDs (folder names) 
% example: celltype = 'VIP'; mice_IDs= {'CL136', 'CL144', 'CL146', 'CL188', 'CL189'};
% returns a struct with high, low and zero reliability ROIs based on 50% of
% CDF for tone and reward responses 
%
% function [p_tone,p_reward] = compare_day1and7_ROI_reliability(celltype,mice_IDs)

clear
close all
%% Enter info for group you want to analyze

celltype = 'PN';


%%
days = {'Day 1', 'Day 7'}; 
outer_folder = 'D:\2P data';
if strcmp(celltype, 'PN')
    c = [0 0 0];
    c1 = c + [0.6 0.6 0.6];
    mice_IDs = {'CL174', 'CL175', 'CL176', 'CL181', 'CL182', 'CL184'};
end
if strcmp(celltype, 'PV')
    mice_IDs = {'CL147','CL172', 'CL195', 'CL196', 'CL198','CL158'};
    c = [0.9 0 0];
    c1 = [1 0.6 0.6];
end
if strcmp(celltype, 'VIP')
    c = [0 0.5 1];
    c1 = [0.6 0.8 1];
    mice_IDs = {'CL136', 'CL144', 'CL146', 'CL188'};
end
if strcmp(celltype, 'SOM')
    c = [0 100/255 0];
    c1 = [0/255 200/255 0/255];
    mice_IDs = {'CL170', 'CL173', 'CL186', 'CL187', 'CL191', 'CL192', 'CL193'};
end

cell_path = strcat(outer_folder, filesep, celltype);

for n = 1:length(mice_IDs)

        %% load data and run response_reliability to get each ROI's response percentage
         current_mouse = mice_IDs{n};
         path = strcat(cell_path, filesep, current_mouse);
%          disp (strcat('Analysing ', current_mouse))
            disp(path)
         [percent_tone,percent_reward] = updated_response_reliability(path,days);
    percent_tone_cell{n} = percent_tone;
    percent_reward_cell{n} = percent_reward;
end

for d = 1:length(days)
    for n = 1:length(mice_IDs)
         num_rois(n) = length(percent_tone);    % number of ROIs for each mouse

        %% pool response reliabilities across mice
        if n == 1 
            percent_tone_pooled_oneday = percent_tone_cell{n}(d,:); 
            percent_reward_pooled_oneday = percent_reward_cell{n}(d,:);
        end
        
        if n > 1
            percent_tone_pooled_oneday = horzcat(percent_tone_pooled_oneday,percent_tone_cell{n}(d,:)); 
            percent_reward_pooled_oneday = horzcat(percent_reward_pooled_oneday,percent_reward_cell{n}(d,:));
        end

    end
    percent_tone_pooled{d} = percent_tone_pooled_oneday;
    percent_reward_pooled{d} = percent_reward_pooled_oneday;
    percent_tone_pooled_oneday = [];
    percent_reward_pooled_oneday = [];
end


        
%% Create CDF
for d = 1:length(days)
    
    temp_tone = percent_tone_pooled{d};
    temp_tone((isnan(temp_tone))) = [];
    percent_tone_pooled_sorted{d} = sort(temp_tone,'ascend');
%     percent_tone_pooled_sorted{d}(percent_tone_pooled_sorted{d}==0) = [];
    
    temp_reward = percent_reward_pooled{d};
    temp_reward((isnan(temp_reward))) = [];
    percent_reward_pooled_sorted{d} = sort(temp_reward,'ascend');
%     percent_reward_pooled_sorted{d}(percent_reward_pooled_sorted{d} ==0) = [];
    
    % removing zeros/ NaNs are unresponsive ROIs could mess with # of ROIs.
    % Safest method is to take a value for each day
    total_rois_pooled_tone(d) = (length(percent_tone_pooled_sorted{d}));
    total_rois_pooled_reward(d) = (length(percent_reward_pooled_sorted{d}));
    
    yt_temp = 1:1:total_rois_pooled_tone(d);
    yt{d} = (yt_temp./total_rois_pooled_tone(d)).*100;
    yr_temp = 1:1:total_rois_pooled_reward(d);
    yr{d} = (yr_temp./total_rois_pooled_reward(d)).*100;
    
    cdf50_tone(d) = median(percent_tone_pooled_sorted{d});
    cdf50_reward(d) = median(percent_reward_pooled_sorted{d});
    
    temp_tone = [];
    temp_reward = [];
end

y_arrow = [50 50]; 
x_arrow_tone = [cdf50_tone(1) cdf50_tone(2)];
x_arrow_reward = [cdf50_reward(1) cdf50_reward(2)];

%% KS test 
[~,p_tone] = kstest2(percent_tone_pooled_sorted{1},percent_tone_pooled_sorted{2});
[~,p_reward] = kstest2(percent_reward_pooled_sorted{1},percent_reward_pooled_sorted{2});
disp(['KS test for tone: p = ', num2str(p_tone)])
disp(['KS test for reward: p = ', num2str(p_reward)])

%% Make figure 
figure(1) 
subplot(1,2,1)  %tone
plot(percent_tone_pooled_sorted{1},yt{1}, 'Color', c1, 'LineWidth', 2);
hold on 
plot(percent_tone_pooled_sorted{2},yt{2}, 'Color', c, 'LineWidth', 2);
title('Tone');
xlabel('Percent of trials an ROI responds to') 
ylabel('Percent of ROIs') 
xlim([0 100])   
SE = [(max(xlim)-8) min(ylim)+8];
text(SE(1), SE(2), ['KS test: p = ', num2str(p_tone)], 'VerticalAlignment','top', 'HorizontalAlignment','right', 'FontSize', 12)
set(gca, 'FontSize', 12)

subplot(1,2,2)  %reward
plot(percent_reward_pooled_sorted{1},yr{1}, 'Color', c1, 'LineWidth', 2);
hold on 
plot(percent_reward_pooled_sorted{2},yr{2}, 'Color', c, 'LineWidth', 2);
title('Reward');
xlabel('Percent of trials an ROI responds to') 
ylabel('Percent of ROIs') 
xlim([0 100])
legend('Day 1', 'Day 7')
SE = [(max(xlim)-8) min(ylim)+8];
text(SE(1), SE(2), ['KS Test: p = ', num2str(p_reward)], 'VerticalAlignment','top', 'HorizontalAlignment','right', 'FontSize', 12)
set(gca, 'FontSize', 12)

suptitle(strcat('Cumulative distribution function - ', celltype))

 