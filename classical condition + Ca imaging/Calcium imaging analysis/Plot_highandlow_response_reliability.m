%updated 2020-09-15 
% This creates line plots of reliabilities for high and low, tone and
% reward

% (1) run updated_classify_day1_ROI_reliability to get the day 1 classifications 
% (2) run updated_response_reliability on day 1 and 7 to get day 1 and 7 reliability values 
% (3) plot the reliability within each class

clear 
close all
clc
%% Enter info for group you want to analyze
% celltype = 'PN';
% mice_IDs ={'CL174', 'CL175', 'CL176', 'CL181', 'CL182', 'CL184'};
% % 
% celltype = 'PV';
% mice_IDs = {'CL147','CL172', 'CL195', 'CL196', 'CL198','CL158'};
% % % 
celltype = 'VIP';
mice_IDs = {'CL136', 'CL144', 'CL146', 'CL188'};

% celltype = 'SOM';
% mice_IDs = {'CL170', 'CL173', 'CL186', 'CL187', 'CL191', 'CL192', 'CL193'};

n = length(mice_IDs);
save = 'N'; % do you want to save this
save_path = 'C:\Users\';
avg_with_mean = 1; % Make this a one if you want to use mean
avg_with_median = 0; % make this a one if you want to use median

%% Get Day 1 classifications
day1 = {'Day 1'};
[day1_reliability] = updated_classify_day1_ROI_reliability(celltype,mice_IDs,day1);
days = {'Day 1', 'Day 7'};

outer_folder = 'D:\2P data';
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
    c1 = [0/255 200/255 0/255];
    c2 = [0 100/255 0];
end

%% Now track those ROIs across days
% day1_reliabiliy holds the indices of the ROIs in different classes
% (classes are high, low and zero reliability based on day 1)

% We will divide up the neurons based on day 1 reliability, and take their
% percent reliability across days

% the new cell arrays will have a cell for each mouse, and within each
% cell, a row is a day, a column is a different ROI but indexing has been
% lost after this
for mouse = 1:n
    current_mouse = mice_IDs{mouse};
    path = strcat(outer_folder, filesep, celltype, filesep, current_mouse);
    disp(path)
    
    %% get reliability values for each class
    percent_tone = [];
    percent_reward = [];
    [percent_tone,percent_reward] = updated_response_reliability(path,days);
    
    %% tone
    % high_rel_tone holds the % reliability over all days for the ROIs that
    % were high reliability on day 1
    curr_index = [];
    curr_index = cell2mat(day1_reliability.tone.high_rel_rois(1, mouse));
    high_rel_tone{mouse} = percent_tone(:,curr_index);
    
    % Get low reliability to tone
    curr_low_index = [];
    curr_low_index = cell2mat(day1_reliability.tone.low_rel_rois(1, mouse));
   
    % get unresponsive to tone
    curr_unres_index = [];
    curr_unres_index = cell2mat(day1_reliability.tone.zero_rel_rois(1, mouse));
    
    % combine the low and unresponsive indices into new low group
    curr_lowzero_index = horzcat(curr_low_index, curr_unres_index);
    low_rel_tone{mouse} = percent_tone(:,curr_lowzero_index);

    %% reward
    %high reliability to reward 
    curr_index = [];
    curr_index = cell2mat(day1_reliability.reward.high_rel_rois(1, mouse));
    high_rel_reward{mouse} = percent_reward(:,curr_index);
          
    % Get low reliability to reward
    curr_low_index = [];
    curr_low_index = cell2mat(day1_reliability.reward.low_rel_rois(1, mouse));
   
    % get unresponsive to reward
    curr_unres_index = [];
    curr_unres_index = cell2mat(day1_reliability.reward.zero_rel_rois(1, mouse));
    
    % combine the low and unresponsive indices into new low group
    curr_lowzero_index = horzcat(curr_low_index, curr_unres_index);
    low_rel_reward{mouse} = percent_reward(:,curr_lowzero_index);
    
    if avg_with_median == 1
        median_rel.tone.high(mouse,:) = nanmedian(high_rel_tone{1,mouse},2);
        median_rel.tone.low(mouse,:) = nanmedian(low_rel_tone{1,mouse},2);
        median_rel.reward.high(mouse,:) = nanmedian(high_rel_reward{1,mouse},2);
        median_rel.reward.low(mouse,:) = nanmedian(low_rel_reward{1,mouse},2);
    end
    if avg_with_mean == 1
        median_rel.tone.high(mouse,:) = nanmean(high_rel_tone{1,mouse},2);
        median_rel.tone.low(mouse,:) = nanmean(low_rel_tone{1,mouse},2);
        median_rel.reward.high(mouse,:) = nanmean(high_rel_reward{1,mouse},2);
        median_rel.reward.low(mouse,:) = nanmean(low_rel_reward{1,mouse},2);
    end        
end

%% Calculate SEM
for d = 1:length(days) 
    
    % in case n changes on different days
    daily_n.tone_high(d) = length(find(~isnan(median_rel.tone.high(:,d))));
    daily_n.tone_low(d) = length(find(~isnan(median_rel.tone.low(:,d))));
    
    daily_n.reward_high(d) = length(find(~isnan(median_rel.reward.high(:,d))));
    daily_n.reward_low(d) = length(find(~isnan(median_rel.reward.low(:,d))));
    
    % Now calculate SEM
    sem.tone_high = nanstd(median_rel.tone.high,0,1)/sqrt(daily_n.tone_high(d));
    sem.tone_low = nanstd(median_rel.tone.low,0,1)/sqrt(daily_n.tone_low(d));
    
    sem.reward_high = nanstd(median_rel.reward.high,0,1)/sqrt(daily_n.reward_high(d));
    sem.reward_low = nanstd(median_rel.reward.low,0,1)/sqrt(daily_n.reward_low(d));

end

%% Calculate mean between mice to get daily mean
overall_avg.tone_high = nanmean(median_rel.tone.high,1);
overall_avg.tone_low = nanmean(median_rel.tone.low,1);

overall_avg.reward_high = nanmean(median_rel.reward.high,1);
overall_avg.reward_low = nanmean(median_rel.reward.low,1);


%% Do t-test to compare day 1 and 7
[~,p.tonehigh] = ttest(median_rel.tone.high(:,1),median_rel.tone.high(:,2)); 
[~,p.tonelow] = ttest(median_rel.tone.low(:,1), median_rel.tone.low(:,2)); 

[~,p.rewardhigh] = ttest(median_rel.reward.high(:,1),median_rel.reward.high(:,2)); 
[~,p.rewardlow] = ttest(median_rel.reward.low(:,1), median_rel.reward.low(:,2)); 

%%
s1 = strcat( celltype, {'-IN Tone High Reliability: Day 1: '}, num2str(round(overall_avg.tone_high(1),2)),{' ± '}, num2str(round(sem.tone_high(1),2)), '%'); 
s2 = strcat( {', Day 7: '}, num2str(round(overall_avg.tone_high(2),2)), {' ± '}, num2str(round(sem.tone_high(2),2)), '%');
disp(strcat(s1,s2));

s1 = strcat( celltype, {'-IN Tone Low Reliability: Day 1: '}, num2str(round(overall_avg.tone_low(1),2)),{' ± '}, num2str(round(sem.tone_low(1),2)), '%'); 
s2 = strcat( {', Day 7: '}, num2str(round(overall_avg.tone_low(2),2)), {' ± '}, num2str(round(sem.tone_low(2),2)), '%');
disp(strcat(s1,s2));

s1 = strcat( celltype, {'-IN Reward High Reliability: Day 1: '}, num2str(round(overall_avg.reward_high(1),2)),{' ± '}, num2str(round(sem.reward_high(1),2)), '%'); 
s2 = strcat( {', Day 7: '}, num2str(round(overall_avg.reward_high(2),2)), {' ± '}, num2str(round(sem.reward_high(2),2)), '%');
disp(strcat(s1,s2)); 

s1 = strcat( celltype, {'-IN Reward Low Reliability: Day 1: '}, num2str(round(overall_avg.reward_low(1),2)),{' ± '}, num2str(round(sem.reward_low(1),2)), '%'); 
s2 = strcat( {', Day 7: '}, num2str(round(overall_avg.reward_low(2),2)), {' ± '}, num2str(round(sem.reward_low(2),2)), '%');
disp(strcat(s1,s2)); 

%% Make tone figure
x = 1:1:length(days);
y_lim_vec = [0 70];
y_star = 60;
close all
figure(1) 

% TONE high reliability ----------------------------------------------
for mouse = 1:n 
    subplot(1,2,1) 
    plot(x,median_rel.tone.high(mouse,:), '-o',...
    'Color', c1,...
    'MarkerFaceColor', c1,...
    'MarkerSize', 3,...
    'MarkerEdgeColor', c1);
    hold on
end
title('High Reliability')
e = errorbar(x, overall_avg.tone_high, sem.tone_high, '-o', 'MarkerSize', 3,...
    'MarkerFaceColor', c2,... 
    'MarkerEdgeColor', c2,...
    'LineWidth', 1.5);
e.Color = c2;
xticks(1:1:length(days))
xticklabels(days)
xlim([0 length(days)+1])
ylabel('Reliability (%)')
ylim(y_lim_vec)
set(gca, 'FontSize', 12)
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['p=',num2str(p.tonehigh)], 'VerticalAlignment','top', 'HorizontalAlignment','right');
if p.tonehigh < 0.05
    text(max(xlim)/2.2, y_star,'*','FontSize',20)
end
%Tone low reliability -----------------------------------------------------
for mouse = 1:n
    subplot(1,2,2)
    plot(x,median_rel.tone.low(mouse,:), '-o',...
    'Color', c1,...
    'MarkerFaceColor', c1,...
    'MarkerSize', 3,...
    'MarkerEdgeColor', c1);
    hold on
end
title('Low Reliability')
e = errorbar(x, overall_avg.tone_low, sem.tone_low, '-o', 'MarkerSize', 3,...
    'MarkerFaceColor', c2,... 
    'MarkerEdgeColor', c2,...
    'LineWidth', 1.5);
e.Color = c2;
xticks(1:1:length(days))
xticklabels(days) 
xlim([0 length(days)+1])
ylim(y_lim_vec)
set(gca, 'FontSize', 12)
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['p = ',num2str(p.tonelow)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
if p.tonelow < 0.05
    text(max(xlim)/2.2, y_star,'*','FontSize',20)
end
suptitle(strcat(celltype, ' - Tone'));

% save
if strcmp(save, 'Y') == 1
    file_name = strcat(celltype,'_tone_reliability_plots.eps');
    s = strcat(save_path, filesep, file_name);
    print(gcf,'-depsc','-painters',s)
    
    file_name = strcat(celltype,'_tone_reliability_plots.tif');
    s = strcat(save_path, filesep, file_name);
    saveas(gcf,s)   % TO SAVE AS TIFF
end
hold off
%% make reward figure
figure(2) 

% high reliability -------------------------------------------
for mouse = 1:n
    subplot(1,2,1) 
    plot(x,median_rel.reward.high(mouse,:), '-o',...
    'Color', c1,...
    'MarkerFaceColor', c1,...
    'MarkerSize', 3,...
    'MarkerEdgeColor', c1);
    hold on
end
title('High Reliability')
e = errorbar(x, overall_avg.reward_high, sem.reward_high, '-o', 'MarkerSize', 3,...
    'MarkerFaceColor', c2,... 
    'MarkerEdgeColor', c2,...
    'LineWidth', 1.5);
e.Color = c2;
xticks(1:1:length(days))
xticklabels(days) 
xlim([0 length(days)+1])
ylabel('Reliability (%)')
ylim(y_lim_vec)
set(gca, 'FontSize', 12)
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['p = ',num2str(p.rewardhigh)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
if p.rewardhigh < 0.05
    text(max(xlim)/2.2, y_star,'*','FontSize',20)
end
% Low reliability -----------------------------------------------
for mouse = 1:n
    subplot(1,2,2)
    plot(x,median_rel.reward.low(mouse,:), '-o',...
    'Color', c1,...
    'MarkerFaceColor', c1,...
    'MarkerSize', 3,...
    'MarkerEdgeColor', c1);
    hold on
end
e = errorbar(x, overall_avg.reward_low, sem.reward_low, '-o', 'MarkerSize', 3,...
    'MarkerFaceColor', c2,... 
    'MarkerEdgeColor', c2,...
    'LineWidth', 1.5);
e.Color = c2;
xticks(1:1:length(days))
xticklabels(days) 
xlim([0 length(days)+1]) 
title('Low reliability')
ylim(y_lim_vec)
set(gca, 'FontSize', 12)
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['p = ', num2str(p.rewardlow)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
set(gca, 'FontSize', 12)
if p.rewardlow < 0.05
    text(max(xlim)/2.2, y_star,'*','FontSize',20)
end
suptitle(strcat(celltype,'- Reward'));

if strcmp(save, 'Y') == 1
    file_name = strcat(celltype,'_reward_reliability_plots.eps');
    s = strcat(save_path, filesep, file_name);
    print(gcf,'-depsc','-painters',s)
    
    file_name = strcat(celltype,'_reward_reliability_plots.tif');
    s = strcat(save_path, filesep, file_name);
    saveas(gcf,s)   % TO SAVE AS TIFF
end

