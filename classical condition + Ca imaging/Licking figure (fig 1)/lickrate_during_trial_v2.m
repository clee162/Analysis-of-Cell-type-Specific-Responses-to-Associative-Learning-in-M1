%plot average lick rate 

clc
clear
close all

post_window = 10; % after tone
pre_window = 2;     % before tone

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

rate_denominator = 0.5;    %in seconds

c1(1,:) = [0.5 0.5 0.5]; %black
c1(2,:) = [1 0.6 0.6]; %light red
c1(3,:) = [153/255 204/255 1];
c1(4,:) = [185/255 236/255 163/255];

c2(1,:) = [0 0 0];
c2(2,:) = [153/255 0 0];%dark red
c2(3,:) = [0 102/255 204/255];  %dark blue
c2(4,:) = [0 102/255 0];    %dark green


for celltype = 1:length(celltypes)
    curr_celltype = celltypes{celltype};
    curr_celltype_mice_IDs = mice_IDs{celltype};
    n(celltype) = length(curr_celltype_mice_IDs);
    
    for d = 1:length(days) 
        mouse_counter = 1;
        curr_day = days{d};
        
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
                    
                    % find licks within our pre to post window
                    lick_ind = find( (lick_start>(tone_start(rew_num)-pre_window)) & (lick_start<(tone_start(rew_num)+post_window)) );
                    if isempty(lick_ind) 
                        lick_time = NaN;
                    else          
                        % relative to tone time 
                        lick_time = lick_start(lick_ind);
                        curr_latencies = lick_time - tone_start(rew_num);
                        counter = -pre_window;
                        for sliding_window = 1:(pre_window+post_window)/rate_denominator
                               curr_licks_ind = find((curr_latencies>counter) & (curr_latencies<(counter+rate_denominator)));
                               if isempty(curr_licks_ind)
                                   curr_licks = 0;
                               else
                                   curr_licks = length(curr_licks_ind);
                               end
                               trial_lick_rate(rew_num,sliding_window) = curr_licks/rate_denominator;
                               counter = counter+rate_denominator;
                        end 
                        trial_mean_anti_lickrate(rew_num) = mean(trial_lick_rate(rew_num,(round(pre_window/rate_denominator)):(round((pre_window+2.5)/rate_denominator))));
                    end 
                end
                lick_rate{d,celltype}(mouse,:) = mean(trial_lick_rate,1);
                all_mice{d}(mouse_counter,:) = mean(trial_lick_rate,1);    % this is for plot with all mice pooled
                mouse_counter = mouse_counter + 1;
                mean_anti_licking{celltype}(mouse,d) = mean(trial_mean_anti_lickrate);
        end
            sem_anti_licking{celltype}(d) = std(mean_anti_licking{celltype}(:,d))/sqrt(n(celltype));
    end
  end
  
  x = (-pre_window+rate_denominator:rate_denominator:post_window);

for d = 1:2 
    mouse_counter = 1;
    for celltype = 1:length(celltypes)
        for mouse = 1:n(celltype)
            all_mice{d}(mouse_counter,:) = lick_rate{d,celltype}(mouse,:);
            all_mice_mean_anti_licking(mouse_counter,d) = mean_anti_licking{celltype}(mouse,d);
            mouse_counter = mouse_counter+1;
        end
    end
end

%% ANOVA for all mice comparing day 1 and 7
[mice,~] = size(all_mice{1});
anticipatory_bins = [(pre_window/rate_denominator):((pre_window+2.5)/rate_denominator)];

for d = 1:2
    for m = 1:mice
        if m == 1 && d == 1
            y = (all_mice{d}(m,anticipatory_bins)');
        else
        y = vertcat(y,all_mice{d}(m,anticipatory_bins)');
        end
    end
end

t = [0:0.5:2.5];
time_grouping = repmat(t,1,mice*2);

d1_string = cell(1,mice*length(t));
d1_string(:) = {'Day 1'};
d7_string = cell(1,mice*length(t));
d7_string(:) = {'Day 7'};
day_grouping = horzcat(d1_string,d7_string);

[p_value,~,stats] = anovan(y,{time_grouping,day_grouping},'varnames',{'time bin','day'},'display', 'off');
% [c,~,~] = multcompare(stats, 'Dimension', [1 2]);

%% figure 
figure(1)
for d = 1:length(days)
    
    sem_all(d,:) = std(all_mice{d})/sqrt(sum(n));
    
    plot(x,mean(all_mice{d}))
    
    if d == 1
           shadedErrorBar_v2(x,mean(all_mice{d}),sem_all(d,:),(c1(1,:)))
       else
            shadedErrorBar_v2(x,mean(all_mice{d}),sem_all(d,:),(c2(1,:)))
       end
       hold on
end

p = patch([0 1 1 0],[0 0 20 20], 'red');
p.FaceColor = [170/255 170/255 170/255];
p.FaceAlpha = 0.2;
p.LineStyle = 'none';

vline(2.5)
xlim([-pre_window post_window])
ylim([0 12])
ylabel('lick rate (per second)')
xlabel('Time from tone onset (s)') 

sem_all_mice_anti_licking = std(all_mice_mean_anti_licking)/sqrt(sum(n));

NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['Day 1 vs 7: p=',num2str(p_value(2))], 'VerticalAlignment','top', 'HorizontalAlignment','right');


save_path = 'C:\Users\clee162\Desktop\Figures\21-05-10\lick plots';
file_name_e = 'Day1and7_meanlickrate.eps';
file_name_t = 'Day1and7_meanlickrate.tif';
s = strcat(save_path, filesep, file_name_e);
print(gcf,'-depsc','-painters',s)  
s = strcat(save_path, filesep, file_name_t);
saveas(gcf,s)   % TO SAVE AS TIFF
