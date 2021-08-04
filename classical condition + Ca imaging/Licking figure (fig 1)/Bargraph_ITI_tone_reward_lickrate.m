%Make a bar graph comparing the lick rate during the first 2.5 s of ITI
%lick bouts with tone lick rate and reward lick rate 

clear
close all
clc

save = 'n';

%% Data path 
window = 2.5; % after tone

PN_mice_IDs ={'CL174', 'CL175', 'CL176', 'CL181', 'CL182', 'CL184'};
PV_mice_IDs = {'CL147','CL172', 'CL195', 'CL196', 'CL198','CL158'}; 
VIP_mice_IDs = {'CL136', 'CL144', 'CL146', 'CL188'};
SOM_mice_IDs = {'CL170', 'CL173', 'CL186', 'CL187', 'CL191', 'CL192', 'CL193'};

outer_folder = 'D:\2P data';
celltypes = {'PN','PV','VIP','SOM'};
mice_IDs{1} = PN_mice_IDs;
mice_IDs{2} = PV_mice_IDs;
mice_IDs{3} = VIP_mice_IDs;
mice_IDs{4} = SOM_mice_IDs;
days = {'Day 1','Day 7'};
t_between_bouts = 3;
fr = 30;

rate_dt = 0.5;    %in seconds

for d = 1:length(days) 
    mouse_counter = 1;
    curr_day = days{d};
        
    for celltype = 1:length(celltypes)
        curr_celltype = celltypes{celltype};
        curr_celltype_mice_IDs = mice_IDs{celltype};
        n(celltype) = length(curr_celltype_mice_IDs);

        for mouse = 1:n(celltype)
            
            ITI_mean_lickrate_trial = [];
            
            %% load data
                curr_mouse = curr_celltype_mice_IDs{mouse};
                path = strcat(outer_folder, filesep, curr_celltype,filesep,curr_mouse,filesep,curr_day);
                disp(path)
                matfiles = dir(fullfile(path,'*mat'));
                nfiles = length(matfiles);
                for j = 1:length(matfiles)
                    s = (fullfile(path,matfiles(j).name));
                    load(s)
                end
                num_trials = length(wt_start);
                [~,frames] = size(df);
                
                for trial = 1:num_trials
                    tone_latencies = [];
                    reward_latencies = [];
                    
                    %% Find ITI lick times
                    ITI_lick_ind = [];
                    lick_bout_start_ind = [];

                    if trial < num_trials    % if it's not the last trial
                        ITI_lick_ind = find((lick_start > (wt_start(trial) + 20)) & (lick_start < (tone_start(trial+1)-window)));
                    end
                    if trial == num_trials    % if it's the last trial
                        ITI_lick_ind = find((lick_start > (wt_start(trial) + 20)) & (lick_start < (frames/fr-window)));
                    end

                    % Find start time of lick bouts
                    curr_trial_ITI_lick_times = lick_start(ITI_lick_ind); %lick times of all ITI lick bouts this trial
                    padded_curr_trial_ITI_lick_times = vertcat(0, curr_trial_ITI_lick_times, 0);
                    diff_lick = diff(padded_curr_trial_ITI_lick_times); % time in between licks
                    lick_bout_start_ind = find(diff_lick > t_between_bouts);  % find where different lick bouts start
                    
                    % find licks within 2.5 s of first lick in bout
                    for j = 1:length(lick_bout_start_ind) 
                        ITI_lickbout_start_time = curr_trial_ITI_lick_times(lick_bout_start_ind(j)); 
                        curr_ITI_lickbout_licks = curr_trial_ITI_lick_times(find( (curr_trial_ITI_lick_times>=ITI_lickbout_start_time) & (curr_trial_ITI_lick_times<(ITI_lickbout_start_time+window)) ));
                        
                        %get times relative to lickbout start
                        curr_ITI_lickbout_licks = (curr_ITI_lickbout_licks) - (ITI_lickbout_start_time);
                        
                        counter = 0;
                        %% calculate lick rate during ITI lickbouts 
                        for sliding_window = 1:window/rate_dt
                            curr_bin_ITI_licks = find((curr_ITI_lickbout_licks>counter) & (curr_ITI_lickbout_licks<(counter+rate_dt)));
                            
                            if isempty(curr_bin_ITI_licks) 
                                curr_ITI_licks = 0;
                            else
                                curr_ITI_licks = length(curr_bin_ITI_licks);
                            end
                            trial_ITI_lick_rate(sliding_window) = curr_ITI_licks/rate_dt;
                            counter = counter+rate_dt;
                        end
                        
                        %  store the mean lickrate for every ITI lick
                        % bout
                        ITI_mean_lickrate_trial{trial}(j) = mean(trial_ITI_lick_rate);    
                    end   
                    
                    %% find licks during tone and reward period 
                    tone_lick_ind = find( (lick_start>(tone_start(trial))) & (lick_start<(tone_start(trial)+window)) );
                    reward_lick_ind = find( (lick_start>(wt_start(trial))) & (lick_start<(wt_start(trial)+window)) );
                    
                    if isempty(tone_lick_ind) 
                        tone_lick_time = NaN;
                    else
                      tone_lick_time = lick_start(tone_lick_ind);
                    end

                    if isempty(reward_lick_ind) 
                        reward_lick_time = NaN;
                    else
                        reward_lick_time = lick_start(reward_lick_ind);
                    end
                    
                    if ~isempty(tone_lick_ind) || ~isempty(reward_lick_ind)                  
                        % relative to onset time  
                        tone_latencies = tone_lick_time - tone_start(trial);
                        reward_latencies = reward_lick_time - wt_start(trial);
                        
                        counter = 0;
                        % calculate lick rate during tone and reward 
                        for sliding_window = 1:window/rate_dt
                            
                               %find the # of licks in each time bin
                               %specified by rate_dt 
                               curr_tone_licks_ind = find((tone_latencies>counter) & (tone_latencies<(counter+rate_dt)));
                               curr_reward_licks_ind = find((reward_latencies>counter) & (reward_latencies<(counter+rate_dt)));

                               if isempty(curr_tone_licks_ind)
                                   curr_tone_licks = 0;                                  
                               else
                                   curr_tone_licks = length(curr_tone_licks_ind);
                               end
                               if isempty(curr_reward_licks_ind)
                                   curr_reward_licks = 0;
                               else
                                   curr_reward_licks = length(curr_reward_licks_ind);
                               end
                               %lick rate during time bin
                               trial_tone_lick_rate(sliding_window) = curr_tone_licks/rate_dt;
                               trial_reward_lick_rate(sliding_window) = curr_reward_licks/rate_dt;
                               counter = counter+rate_dt;
                        end 
                        
                        tone_mean_lickrate_trial(trial) = mean(trial_tone_lick_rate);
                        reward_mean_lickrate_trial(trial) = mean(trial_reward_lick_rate);
                        
                        end 
                   
                end
                %% HERE ARE YOUR MEAN VALUES              
                all_ITI_lickrates = cell2mat(ITI_mean_lickrate_trial);
                all_mice_mean_ITI(mouse_counter,d) = mean(all_ITI_lickrates);
                all_mice_mean_tone(mouse_counter,d) = mean(tone_mean_lickrate_trial);    % this is for plot with all mice pooled
                all_mice_mean_reward(mouse_counter,d) = mean(reward_mean_lickrate_trial);
                mouse_counter = mouse_counter + 1;
                

    end
    end
                 
end
%% SEM
sem_tone_lickrate= std(all_mice_mean_tone)/sqrt(sum(n));
sem_reward_lickrate = std(all_mice_mean_reward)/sqrt(sum(n));
sem_ITI_lickrate = std(all_mice_mean_ITI)/sqrt(sum(n));

%% make figure 
a = [0 0.12 -0.1 -0.02 0.04 -0.09 0.12 -0.05 0.03];
r = -0.12 + (0.12+0.12)*rand(mouse_counter,1);
x = 1:1:3;

%% arrange data for easier graphing 
d1_bardata = horzcat((all_mice_mean_ITI(:,1)),(all_mice_mean_tone(:,1)), (all_mice_mean_reward(:,1)));
d2_bardata = horzcat((all_mice_mean_ITI(:,2)),(all_mice_mean_tone(:,2)), (all_mice_mean_reward(:,2)));
d1_sem = horzcat(sem_ITI_lickrate(1), sem_tone_lickrate(1), sem_reward_lickrate(1));
d2_sem = horzcat(sem_ITI_lickrate(2), sem_tone_lickrate(2), sem_reward_lickrate(2));

%% ANOVA 
[p1,~,stats] = anova1(d1_bardata, [], 'off');
[c1,~,~] = multcompare(stats);
[p2,~,stats] = anova1(d2_bardata, [], 'off');
[c2,~,~] = multcompare(stats,'display', 'off');
close all

%% display mean and sem
disp(strcat({' Day 1 ITI: '}, num2str(mean(d1_bardata(:,1))), {' ± '}, num2str(d1_sem(1))));
disp(strcat({' Day 1 Tone: '}, num2str(mean(d1_bardata(:,2))), {' ± '}, num2str(d1_sem(2))));
disp(strcat({' Day 1 Reward: '}, num2str(mean(d1_bardata(:,3))), {' ± '}, num2str(d1_sem(3))));

disp(strcat({' Day 7 ITI: '}, num2str(mean(d2_bardata(:,1))), {' ± '}, num2str(d2_sem(1))));
disp(strcat({' Day 7 Tone: '}, num2str(mean(d2_bardata(:,2))), {' ± '}, num2str(d2_sem(2))));
disp(strcat({' Day 7 Reward: '}, num2str(mean(d2_bardata(:,3))), {' ± '}, num2str(d2_sem(3))));

%% Day 1 figure
figure(1)
bar(x,mean(d1_bardata),'FaceColor',[1 1 1])
hold on 
er = errorbar(x,mean(d1_bardata),d1_sem,d1_sem);
er.Color = [0 0 0];
er.LineStyle = 'none';  
hold on
sz = 25;

for i = 1:3
    for point = 1:mouse_counter-1
        scatter_x = r(point)+i;
        scatter(scatter_x,d1_bardata(point,i),sz,...
            'MarkerFaceAlpha', 0.8,...
            'MarkerEdgeColor', 'k')
        hold on
    end
end
xlim([0 4])
% ylim([0 50])
ylabel('Lick Rate (licks/second)')
title('Day 1 Lick Rate')
xticks(1:3)
xticklabels({'ITI', 'Tone', 'Reward'})
NE = [1.5 max(ylim)*0.95];
text(NE(1), NE(2), [{'ITI vs tone p=', num2str(c1(1,6))}, {'ITI vs reward p=' num2str(c1(2,6))}, {'Tone vs reward p =', num2str(c1(3,6))}], 'VerticalAlignment','top', 'HorizontalAlignment','right');

figure(2)
bar(x,mean(d2_bardata),'FaceColor',[1 1 1])
hold on 
er = errorbar(x,mean(d2_bardata),d2_sem,d2_sem);
er.Color = [0 0 0];
er.LineStyle = 'none';  
hold on
sz = 25;

for i = 1:3
    for point = 1:mouse_counter-1
        scatter_x = r(point)+i;
        scatter(scatter_x,d2_bardata(point,i),sz,...
            'MarkerFaceAlpha', 0.8,...
            'MarkerEdgeColor', 'k')
        hold on
    end
end
xlim([0 4])
% ylim([0 50])
ylabel('Lick Rate (licks/second)')
title('Day 7 Lick Rate')
xticks(1:3)
xticklabels({'ITI', 'Tone', 'Reward'})
hold on
figure(2)
% NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
NE = [1.5 max(ylim)*0.95];
text(NE(1), NE(2), [{'ITI vs tone p=', num2str(c2(1,6))}, {'ITI vs reward p=' num2str(c2(2,6))}, {'Tone vs reward p =', num2str(c2(3,6))}], 'VerticalAlignment','top', 'HorizontalAlignment','right');


if strcmp(save, 'Y') == 1
    figure(1) 
    save_path = 'C:\Users\clee162\Desktop\Figures\21-05-10\lick plots';
    file_name_t = 'Bargraph_ITI_Tone_Reward_mean_lickrate_D1.tif';
    file_name_e = 'Bargraph_ITI_Tone_Reward_mean_lickrate_D1.eps';
    s = strcat(save_path, filesep, file_name_e);
    print(gcf,'-depsc','-painters',s)  
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
    % 
    figure(2) 
    save_path = 'C:\Users\clee162\Desktop\Figures\21-05-10\lick plots';
    file_name_e = 'Bargraph_ITI_Tone_Reward_mean_lickrate_D7.eps';
    file_name_t = 'Bargraph_ITI_Tone_Reward_mean_lickrate_D7.tif';
    s = strcat(save_path, filesep, file_name_e);
    print(gcf,'-depsc','-painters',s) 
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
end