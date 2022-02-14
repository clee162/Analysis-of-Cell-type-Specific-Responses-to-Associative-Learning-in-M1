% make a plot showing lick bouts during the ITI

clear
close all
clc

%% inputs
celltype = 'PN';
mouse = 'CL176';
day = 'Day 1';
save = 'l';
outer_folder = 'D:\2P data\';
save_path = 'C:\Users\'; 
save = 'n';

%% if you want to specify trials 
trials_for_plot = 1:11;
num_trials_plot = length(trials_for_plot);

%% parameters
window = 2.5; 
t1_plot = 20;       % time in seconds before tone
t2_plot = 5;    %time in s after tone
tone_length = 1; 
delay_length = 1.5; 
% num_trials_plot = 10;
t_between_bouts = 3;   % how much time between licks before it is a new bout

%% load data
path = strcat(outer_folder,celltype,filesep,mouse,filesep,day);
matfiles =  dir(fullfile(path, '*.mat'));
for i = 1:length(matfiles)
    s = (fullfile(path, matfiles(i).name));
    load(s)
end    
    
%% pre-processing 
[~, wt_start, lick_start, tone_start] = daily_preprocessing(df, wt_start, lick_start, tone_start);
num_trials = length(wt_start);

for j = 1:num_trials_plot

    trial = trials_for_plot(j);
   %% find all licks within the span of the figure (from -t1:+t2) - this includes ITI lick bouts and other licks
     trial_lick_ind = find( (lick_start > (tone_start(trial+1)-t1_plot)) & (lick_start < (tone_start(trial+1)+t2_plot)) );

   %% get all lick times relative to the next tone start - this includes ITI lick bouts and other licks. this is what we need to plot 
    trial_lick_times = lick_start(trial_lick_ind);    %lick times of all beam breaks within ITI window  
    trial_lick_times_rel =  trial_lick_times - tone_start(trial+1);
    
    %% find lick bout start and stop
    ITI_lick_times_rel = trial_lick_times_rel(trial_lick_times_rel < -2.5); % find licks that end at least 2.5 s before tone (remember tone = time 0)    
    num_licks_in_ITI = length(ITI_lick_times_rel);   % number of licks in this ITI  
    padded_ITI_lick_times = vertcat(0, ITI_lick_times_rel);    % lick times padded with 0s so we can take diff
    diff_lick = diff(padded_ITI_lick_times);   % this will tell us how far apart lick bouts are
    lickbout_start_ind = find(diff_lick < -t_between_bouts | diff_lick > t_between_bouts);
    if ~isempty(lickbout_start_ind) 
        lickbout_start_times = ITI_lick_times_rel(lickbout_start_ind);
        temp_lickbout_start_ind = lickbout_start_ind;
        temp_lickbout_start_ind(1) = [];
        lickbout_end_ind = temp_lickbout_start_ind - 1;
        lickbout_end_ind = vertcat(lickbout_end_ind, num_licks_in_ITI);
        lickbout_end_times = ITI_lick_times_rel(lickbout_end_ind);
        num_lick_bouts = length(lickbout_start_ind);
    end
    
    %% make plot
    y = ones(1,length(trial_lick_times_rel)).*j; 
    tone_start_plot = 0;
    tone_end_plot = tone_length;
    
    figure(1)
    s = scatter(trial_lick_times_rel, y);
    s.MarkerEdgeColor = 'k';
    % lick bout patches 
    f = [1 2 3 4];
    if ~isempty(lickbout_start_ind) 
        for bout = 1:num_lick_bouts
            v = [lickbout_start_times(bout) (j-0.25); lickbout_end_times(bout) (j-0.25); lickbout_end_times(bout) (j+0.25); lickbout_start_times(bout) (j+0.25)];
            p_l = patch('Faces',f,'Vertices',v,'Facecolor', [0.5 .3 0.6]);
            p_l.FaceAlpha = 0.3;
            p_l.EdgeColor = 'none';
        end
    end
    hold on

    %% with lines instead 
    figure(2)
    a = repelem(trial_lick_times_rel,2);
    y_l = [(j-0.25) (j+0.25)];
    b = repmat(y_l,1,length(trial_lick_times_rel));
    counter = 1;
    for i = 1:length(trial_lick_times_rel)
        line(a(counter:counter+1), b(counter:counter+1), 'Color', [0 0 0])
        counter = counter + 2;
        hold on
    end
    
    % lick bout patches 
    if ~isempty(lickbout_start_ind) 
        for bout = 1:num_lick_bouts
            v = [lickbout_start_times(bout) (j-0.25); lickbout_end_times(bout) (j-0.25); lickbout_end_times(bout) (j+0.25); lickbout_start_times(bout) (j+0.25)];
            p_l = patch('Faces',f,'Vertices',v,'Facecolor', [0.5 .3 0.6]);
            p_l.FaceAlpha = 0.3;
            p_l.EdgeColor = 'none';
        end
    end    
end

% let's come back to the figures outside of the loop
figure(1) 
% tone patch 
v = [tone_start_plot (num_trials_plot+1); tone_end_plot (num_trials_plot+1); tone_end_plot (num_trials_plot+1.5); tone_start_plot (num_trials_plot+1.5)];
p_t = patch('Faces',f,'Vertices',v,'FaceColor', [0.5 0.5 0.5]); 
p_t.EdgeColor = 'none';
axis([-t1_plot t2_plot 0 num_trials_plot+2])
ylabel('Trial')
xlabel('Time from tone onset')
yticks(1:num_trials_plot)
vline(tone_length + delay_length)
title(strcat(mouse, {' '}, day))

figure(2)
v = [tone_start_plot (num_trials_plot+1); tone_end_plot (num_trials_plot+1); tone_end_plot (num_trials_plot+1.5); tone_start_plot (num_trials_plot+1.5)];
p_t = patch('Faces',f,'Vertices',v,'FaceColor', [0.5 0.5 0.5]); 
p_t.EdgeColor = 'none';
axis([-t1_plot t2_plot 0 num_trials_plot+2])
ylabel('Trial')
xlabel('Time from tone onset')
yticks(1:num_trials_plot)
vline(tone_length + delay_length)
title(strcat(mouse, {' '}, day))

if strcmp(save, 'Y') == 1
    figure(1)
    file_name_e = strcat(celltype,'ITI_licktimes_circle.eps'); % .eps OR .tif
    file_name_t = strcat(celltype,'_ITI_licktimes_circle.tif'); % .eps OR .tif
    s = strcat(save_path, filesep, file_name_e);
    print(gcf,'-depsc','-painters',s) % TO SAVE AS EPS
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
    hold off 
    figure(2)
    file_name_e = strcat(celltype,'ITI_licktimes_line.eps'); % .eps OR .tif
    file_name_t = strcat(celltype,'_ITI_licktimes_line.tif'); % .eps OR .tif
    s = strcat(save_path, filesep, file_name_e);
    print(gcf,'-depsc','-painters',s) % TO SAVE AS EPS
    s = strcat(save_path, filesep, file_name_t);
    saveas(gcf,s)   % TO SAVE AS TIFF
end
