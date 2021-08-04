clear 
clc
close all
%clf

celltype = 'PV';
mouse = 'CL172';
day = 'Day 7';
reward_num = 14; 
save = 'n';
save_path = 'C:\Users\';
color_bar = [0 7];
catch_trial = 0;
fr = 30;
dt = 1/fr;
t1 = -1;
t2 = 6;
tone_dur = 1;
delay_length = 1.5;
fig_length = -t1 + t2;

%% Get number of ROIs on day 1
day1_path = strcat('D:\2P data\',celltype,filesep,mouse,filesep,'Day 1');
matfiles =  dir(fullfile(day1_path, '*.mat'));
nfiles = length(matfiles);
for i = 1:length(matfiles)
    s = (fullfile(day1_path, matfiles(i).name));
    load(s)
end  
    
[num_rois_day1, ~ ] = size(df);
df = [];
real_rois = [];
wt_start = [];
tone_start = [];
lick_start = [];

%% load data
path = strcat('D:\2P data\',celltype,filesep,mouse,filesep,day);
matfiles =  dir(fullfile(path, '*.mat'));
    nfiles = length(matfiles);
    for i = 1:length(matfiles)
        s = (fullfile(path, matfiles(i).name));
        load(s)
    end    
    
real_rois(real_rois > num_rois_day1) = [];
real_rois(real_rois == 0) = [];
df = df(real_rois,:);

% Restrict t to time of interest for sorting ( reward -2 and +5 seconds)
%% find frames for the trial
if catch_trial == 0 
    tone_t = tone_start(reward_num);
    fig_start = tone_t+t1; % start time for figure
    fig_end = tone_t + t2; % end time for figure
else
    tone_t = catch_t(catch_trial,2);
    fig_start = tone_t+t1; % start time for figure
    fig_end = fig_start + t2; 
end

ind_s = round(fig_start/dt);      %start index    
ind_e = round(fig_end/dt);      %end index
t_window = fig_start:dt:fig_end;

%% find lick times
[lick_ind,~] = find(lick_start > fig_start & lick_start < fig_end); %find licks during time window
if ~isempty(lick_ind)
    lick_times = lick_start(lick_ind);
    lick_frames_rel = (lick_times - fig_start)*fr;  % lick frame RELATIVE to start of figure for current trial 
end
DF = df(:,ind_s:ind_e);     % raw dF/F for the image period only
[r,n_frames] = size(DF);

%% Do preprocessing and get z-score
[df_z, ~, ~, ~] = daily_preprocessing(df, wt_start, lick_start, tone_start);
DF_z = df_z(:,ind_s:ind_e);

%% sort by max 
DF_z_after_tone = DF_z(:,abs(t1)*fr:n_frames);
DF_after_tone = DF(:,abs(t1)*fr:n_frames);
index_z = zeros(r,2);
index = zeros(r,2);
for i = 1:r
        index_z(i,:) = [i,find(DF_z_after_tone == max(DF_z_after_tone(i,:)))];    %find index of max value following the tone
        index(i,:) = [i,find(DF_after_tone == max(DF_after_tone(i,:)))];    %find index of max value following the tone
end
[a,b]=sort(index_z(:,2));
index_z_sort = [index_z(b,1) a];
DF_z_sorted = DF_z(index_z_sort(:,1),:);

[a,b] = sort(index(:,2));
index_sort = [index(b,1) a];
DF_sorted = DF(index_sort(:,1),:);
%% Create figure

% DF 
figure(1) 
% step (1) - add filler rows
space_h = round(r*.1);                               % number of rows
space_holder = ones(space_h, n_frames)*20;
DF_sorted = vertcat(space_holder, DF_sorted); 
% step (2) - make heatmap
imagesc(DF_sorted);
colormap(hot) 
caxis([0 1])
% step (3) - white patch 
p1 = patch([0 n_frames  n_frames 0], [0 0 space_h space_h] , 'white');
p1.EdgeColor = [1 1 1];
% step (4) - tone patch
p_t = patch([-(t1)*fr (-t1 + tone_dur)*fr (-t1 + tone_dur)*fr -t1*fr],[0 0 round(space_h.*.9) round(space_h*.9)], 'red');
p_t.FaceColor = [100/255 100/255 100/255];
p_t.FaceAlpha = 0.5;
p_t.EdgeColor = [170/255 170/255 170/255];
% step (5) - lick lines
if ~isempty(lick_ind)
    a = repelem(lick_frames_rel,2);             % x values
    y = [0 round(space_h*.9)];               
    b = repmat(y,1,length(lick_frames_rel));    % y values
    counter = 1;
    for i = 1:length(lick_ind)
        line(a(counter:counter+1), b(counter:counter+1), 'Color', [0 0 0])
        counter = counter + 2;
    end
end    
% step (7) - reward line 
vline((-t1 + tone_dur + delay_length)*fr,'w')
title('DF - sorted')
% colorbar
% % set(gca,'Visible','off')
% box off
% set(findall(gca, 'type', 'text'), 'visible', 'on')
% set(gcf, 'InvertHardCopy', 'off');
% save_path = strcat('C:\Users\clee162\Desktop\Figures\20-09-23\DF 0 to 1\', celltype);
% s = strcat(save_path,filesep,celltype,'_',mouse,'_',day,'_',num2str(reward_num),'SORTED');
% if strcmp(save, 'Y') == 1
%     saveas(gcf,s,'tif')
% end


%% DF - z-scored
figure(2)
DF_z_sorted = vertcat(space_holder,DF_z_sorted);
imagesc(DF_z_sorted);
colormap(hot) 
caxis(color_bar) 
% step (3) - white patch 
p1 = patch([0 n_frames n_frames 0], [0 0 space_h+1 space_h+1] , 'white');
p1.EdgeColor = [1 1 1];
% step (4) - tone patch
p_t = patch([-(t1)*fr (-t1 + tone_dur)*fr (-t1 + tone_dur)*fr -t1*fr],[0 0 round(space_h.*.9) round(space_h*.9)], 'red');
p_t.FaceColor = [100/255 100/255 100/255];
p_t.FaceAlpha = 0.5;
p_t.EdgeColor = [170/255 170/255 170/255];
% step (5) - lick lines
if ~isempty(lick_ind)
    counter = 1;
    for i = 1:length(lick_ind)
        line(a(counter:counter+1), b(counter:counter+1), 'Color', [0 0 0])
        counter = counter + 2;
    end
end  
vline((-t1 + tone_dur + delay_length)*fr,'w')
colorbar
title(strcat(celltype,{' '},mouse,{' '},day,{' '},{'trial '},num2str(reward_num)))
% set(gca,'Visible','off')
box off
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gcf, 'InvertHardCopy', 'off');
if strcmp(save, 'Y') == 1
   s = strcat(save_path,filesep,celltype,'_',mouse,'_',day,'_',num2str(reward_num),'_heatmap.eps');
   print(gcf,'-depsc','-painters',s) % TO SAVE AS EPS
   s = strcat(save_path,filesep,celltype,'_',mouse,'_',day,'_',num2str(reward_num),'_heatmap.tif');
   saveas(gcf,s)
end
