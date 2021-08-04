clear
close all
clc

addpath(genpath(pwd))
save = 'n'; % do you want to save this

% PN
% celltype = 'PN';
% path{1} = 'D:\Data - cleaned\PN\CL174';
% path{2} = 'D:\Data - cleaned\PN\CL175';
% path{3} = 'D:\Data - cleaned\PN\CL176';
% path{4} = 'D:\Data - cleaned\PN\CL181';
% path{5} = 'D:\Data - cleaned\PN\CL182';
% path{6} = 'D:\Data - cleaned\PN\CL184';

%% VIP
celltype = 'VIP';
path{1} = 'D:\Data - cleaned\VIP\CL136';
path{2} = 'D:\Data - cleaned\VIP\CL144';
path{3} = 'D:\Data - cleaned\VIP\CL146';
path{4} = 'D:\Data - cleaned\VIP\CL188';

%% PV 
% celltype = 'PV';
% path{1} = 'D:\Data - cleaned\PV\CL147';
% path{2} = 'D:\Data - cleaned\PV\CL172';
% path{3} = 'D:\Data - cleaned\PV\CL195';
% path{4} = 'D:\Data - cleaned\PV\CL196';
% path{5} = 'D:\Data - cleaned\PV\CL198';
% path{6} = 'D:\Data - cleaned\PV\CL158';

% % SOM 
% celltype = 'SOM';
% path{1} = 'D:\Data - cleaned\SOM\CL170';
% path{2} = 'D:\Data - cleaned\SOM\CL173';
% path{3} = 'D:\Data - cleaned\SOM\CL186';
% path{4} = 'D:\Data - cleaned\SOM\CL187';
% path{5} = 'D:\Data - cleaned\SOM\CL191';
% path{6} = 'D:\Data - cleaned\SOM\CL192';
% path{7} = 'D:\Data - cleaned\SOM\CL193';

n = length(path);
fr = 30;
window_l = 2.5;       % in seconds
all_days = {'Day 1','Day 7'};

for d = 1:length(all_days)
    day = all_days{d};
    disp(strcat('Analyzing ', day));  
   
    for mouse = 1:length(path)
        disp(['Analysing mouse ' num2str(mouse)])
         
        df = [];
        tone_start = [];
        wt_start = [];
        real_rois = [];

        %% load data
        curr_path{mouse} = strcat(path{mouse}, '\', day);
        matfiles =  dir(fullfile(curr_path{mouse}, '*.mat'));
        nfiles = length(matfiles);
            for j = 1:length(matfiles)
                s = (fullfile(curr_path{mouse}, matfiles(j).name));
                load(s)
            end          
            
        %% exclude ROIs that don't fire at least once in the entire session
        if d == 1
            [num_d1_rois(mouse),~] = size(df);
        end
       
        [~, frames] = size(df);
        if isempty(find(real_rois))
            disp([strcat(num2str(day) ,' mouse ', num2str(mouse), ' has no responsive rois')])
            daily_avg.tone(mouse,d) = NaN; 
            daily_avg.reward(mouse,d) = NaN;  
        else
            
            unresponsive_ROIs = find(real_rois == 0);
            df = df(1:num_d1_rois(mouse),:);
            if ~isempty(unresponsive_ROIs)
                unresponsive_ROIs(unresponsive_ROIs >= num_d1_rois(mouse)) = [];
                df(unresponsive_ROIs,:) = [];
            end
            [num_rois,~] = size(df);
            
            %% remove tone_starts that are past the imaging
            tone_start = tone_start(tone_start < (frames/fr - (window_l)*fr));
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
            
            window = 2.5;
            [~, median_percent_active_TONE] = find_indices_of_tone_active_ROIs_v3(df_z, tone_start, window);  
            window = 2.5;
            [~, median_percent_active_REWARD] = find_indices_of_reward_active_ROIs_v4(window, df_z, wt_start);
            
            daily_avg.tone(mouse,d) = median_percent_active_TONE; 
            daily_avg.reward(mouse,d) = median_percent_active_REWARD;
        end
    end
end

%% average across all mice 
overall_avg.tone = nanmean(daily_avg.tone,1);
overall_avg.reward = nanmean(daily_avg.reward,1);

%% SEM
for d = 1:length(all_days)
    
    n_tone = length(daily_avg.tone(~isnan(daily_avg.tone(:,d))));
    sem.tone(d) = nanstd(daily_avg.tone(:,d))/sqrt(n_tone);
    
    n_reward = length(daily_avg.reward(~isnan(daily_avg.reward(:,d))));
    sem.reward(d) = nanstd(daily_avg.reward(:,d))/sqrt(n_reward);
end

%% t-test
[~,p.tone] = ttest(daily_avg.tone(:,1),daily_avg.tone(:,2));
[~,p.reward] = ttest(daily_avg.reward(:,1),daily_avg.reward(:,2));

%% make figure

if strcmp(celltype, 'PN')
    c = [0 0 0];
    c1 = c + [0.6 0.6 0.6];
    mice_IDs = {'CL174', 'CL175', 'CL176', 'CL181', 'CL182', 'CL184'};
end
if strcmp(celltype, 'PV')
    mice_IDs = {'CL147', 'CL172', 'CL195', 'CL196'};
    c = [0.9 0 0];
    c1 = [1 0.6 0.6];
end
if strcmp(celltype, 'VIP')
    c = [0 0.5 1];
    c1 = [0.6 0.8 1];
    mice_IDs = {'CL136', 'CL144', 'CL146', 'CL188', 'CL199'};
end
if strcmp(celltype, 'SOM')
    c1 = [0/255 200/255 0/255];
    c = [0 100/255 0];
    mice_IDs = {'CL170', 'CL173', 'CL186', 'CL187', 'CL191', 'CL192', 'CL193'};
end


x = 1:length(all_days);
y_axis_vec = [0 60];
y_star = 40;

s1 = strcat( celltype, {' Tone: Day 1: '}, num2str(round(overall_avg.tone(1),2)),{' ± '}, num2str(round(sem.tone(1),2)), '%'); 
s2 = strcat( {' Day 7 : '}, num2str(round(overall_avg.tone(2),2)), {' ± '}, num2str(round(sem.tone(2),2)), '%');
disp(strcat(s1,s2));

s1 = strcat( celltype, {' Reward: Day 1: '}, num2str(round(overall_avg.reward(1),2)),{' ± '}, num2str(round(sem.reward(1),2)), '%'); 
s2 = strcat( {' Day 7 : '}, num2str(round(overall_avg.reward(2),2)), {' ± '}, num2str(round(sem.reward(2),2)), '%');
disp(strcat(s1,s2)); 
%% 
figure(1) 

% subplot(1,2,1)
for i = 1:length(path) 
    plot(x,daily_avg.tone(i,:), '-o',...
    'Color', c1,... 
   'MarkerFaceColor', c1,...
   'MarkerSize', 4,...
   'MarkerEdgeColor', c1);
    hold on
end
e = errorbar(x, overall_avg.tone, sem.tone, '-o', 'MarkerSize', 4,...
    'MarkerFaceColor', c,... 
    'MarkerEdgeColor', c,...
    'LineWidth', 1.5);
e.Color = c;
box off
xticks(1:1:length(all_days))
xlim([0 length(all_days)+1])
ylim(y_axis_vec) 
xticklabels(all_days)
title(strcat(celltype,'_ Tone - 2.5s window')) 
ylabel('Percent ROIs') 
% xlabel('Day')
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['p=',num2str(p.tone)], 'VerticalAlignment','top', 'HorizontalAlignment','right');
if p.tone < 0.05
    text(max(xlim)/2.2, y_star,'*','FontSize',20)
end
%% save
if strcmp(save, 'Y') == 1
    save_path = 'C:\Users\clee162\Desktop\Figures\21-05-10';
    % EPS vs TIF
    file_name_t = strcat(celltype,'_percent_responsive_TONE.tif'); % .eps OR .tif
    file_name_e = strcat(celltype,'_percent_responsive_TONE.eps');
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
    s = strcat(save_path,filesep, file_name_e);
    print(gcf,'-depsc','-painters',s)     % TO SAVE AS EPS
end
hold off
%% 
figure(2)
% subplot(1,2,2)
for i = 1:length(path) 
    plot(x,daily_avg.reward(i,:), '-o',...
    'Color', c1,... 
   'MarkerFaceColor', c1,...
   'MarkerSize', 4,...
   'MarkerEdgeColor', c1);
    hold on
end
e = errorbar(x, overall_avg.reward, sem.reward, '-o', 'MarkerSize', 4,...
    'MarkerFaceColor', c,... 
    'MarkerEdgeColor', c,...
    'LineWidth', 1.5);
e.Color = c;
box off
xticks(1:1:length(all_days))
xlim([0 length(all_days)+1])
xticklabels(all_days)
title(strcat(celltype,'_ Reward - 2.5s window')) 
ylim(y_axis_vec)           
ylabel('Percent ROIs') 
% xlabel('Day')
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['p=',num2str(p.reward)], 'VerticalAlignment','top', 'HorizontalAlignment','right');
if p.reward < 0.05
   text(max(xlim)/2.2, y_star,'*','FontSize',20)
end

%% save
if strcmp(save, 'Y') == 1
    save_path = 'C:\Users\clee162\Desktop\Figures\21-05-10';
    file_name_t = strcat(celltype,'_percent_responsive_REWARD.tif'); % .eps OR .tif
    file_name_e = strcat(celltype,'_percent_responsive_REWARD.eps');
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
    s = strcat(save_path,filesep, file_name_e);
    print(gcf,'-depsc','-painters',s)     % TO SAVE AS EPS
end

