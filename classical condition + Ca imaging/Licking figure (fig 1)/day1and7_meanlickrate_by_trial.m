%for day 1, plot the mean lick rate during the anticipatory period for all
%trials in day 1

clc
clear
close all

save = 'v';

post_window = 2.5; % after tone

c = {[.6 .6 .6], [0 0 0]};

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

mouse_counter = 1;
rate_denominator = 0.5;    %in seconds
x = 1:33; %x axis of plot (trials)

for d = 1:2
    curr_day = days{d};
    mouse_counter = 1;
    trial_mean_anti_lickrate = [];
    
    for celltype = 1:length(celltypes)
        curr_celltype = celltypes{celltype};
        curr_celltype_mice_IDs = mice_IDs{celltype};
        n(celltype) = length(curr_celltype_mice_IDs);

        for mouse = 1:n(celltype)
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

            for rew_num = 1:length(wt_start)
                    curr_latencies = [];
                    trial_lick_rate = [];
                    % find licks within our pre to post window
                    lick_ind = find( (lick_start>(tone_start(rew_num))) & (lick_start<(tone_start(rew_num)+post_window)) );
                    
                    if isempty(lick_ind) 
                        lick_time = NaN;
                    else 
                    % relative to tone time 
                        lick_time = lick_start(lick_ind);
                        curr_latencies = lick_time - tone_start(rew_num);
                        counter = 0;
                        
                        for sliding_window = 1:post_window/rate_denominator
                               curr_licks_ind = find((curr_latencies>counter) & (curr_latencies<(counter+rate_denominator)));
                               if isempty(curr_licks_ind)
                                   curr_licks = 0;
                               else
                                   curr_licks = length(curr_licks_ind);
                               end
                               trial_lick_rate(sliding_window) = curr_licks/rate_denominator;
                               counter = counter+rate_denominator;
                        end 
                        trial_mean_anti_lickrate{rew_num}(mouse_counter) = mean(trial_lick_rate);                        
                    end
                        

            end
    %         mean_anti_licking{celltype}(mouse,rew_num) = mean(trial_mean_anti_lickrate);
            mouse_counter = mouse_counter + 1;
        end
    end
    
    for rew_num = 1:33
        a = cell2mat(trial_mean_anti_lickrate(rew_num));
        mean_anti_lickrate(rew_num) = mean(a); 
        sem_anti_licking(rew_num) = std(a)/sqrt(length(a));
        
        % prepare data for 2way ANOVA
        cur_trial_grouping = repmat(rew_num,1,length(a)); 
        if d == 1 
            cur_day_grouping = repmat(1,1,length(a));
            if rew_num == 1 
                trial_grouping = cur_trial_grouping;
                day_grouping = cur_day_grouping;
                anova_data = a;
            end
        else 
            cur_day_grouping = repmat(7,1,length(a));
        end
        trial_grouping = horzcat(trial_grouping,cur_trial_grouping);
        day_grouping = horzcat(day_grouping,cur_day_grouping);
        anova_data = horzcat(anova_data,a);
        
    end

    
    figure(1) 
    plot(x,mean_anti_lickrate)
    shadedErrorBar_v2(x,mean_anti_lickrate, sem_anti_licking,c{d})
    hold on
end
xlim([1 33])
% ylim([0 3])
ylabel('lick rate (per second)')
xlabel('Trial number') 
title('Day 1 vs Day 7 lickrate')

[p,~,stats] = anovan(anova_data,{trial_grouping,day_grouping},'varnames',{'trial','day'},'display', 'off');
% [c,~,~] = multcompare(stats, 'Dimension', [1 2]);
figure(1) 
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['2-way ANOVA: trial p =',num2str(p(1)), 'Day p - ', num2str(p(2))], 'VerticalAlignment','top', 'HorizontalAlignment','right');

if strcmp(save, 'Y') == 1
    save_path = 'C:\Users\clee162\Desktop\Figures\21-05-10\lick plots';
    file_name_e = 'Day1and7_meanlickrate_per_trial.eps';
    file_name_t = 'Day1and7_meanlickrate_per_trial.tif';
    s = strcat(save_path, filesep, file_name_e);
    print(gcf,'-depsc','-painters',s)  
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
end
    