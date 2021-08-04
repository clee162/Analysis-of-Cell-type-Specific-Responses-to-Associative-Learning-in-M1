% plots a heat map for an ITI and a reawrd

clear 
clc
close all

%% inputs
celltype = 'PN';
mouse = 'CL175';
day = 'Day 7';
save = 'Y';
trial_plot = 2;
tone_dur = 1;
%% parameters
fr = 30;
dt = 1/fr;
t1 = -1;
t2 = 6;
fig_length = abs(-t1) + t2;
win = 2.5;
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
for i = 1:length(matfiles)
    s = (fullfile(path, matfiles(i).name));
    load(s)
end    
    
real_rois(real_rois > num_rois_day1) = [];
real_rois(real_rois == 0) = [];
df = df(real_rois,:);

%% FIND ITI LICKS

for trial = trial_plot
    ITI_lick_ind = find((lick_start > (wt_start(trial) + 20)) & (lick_start < (tone_start(trial+1)-win)));

    % Find start of lick bouts
    curr_ITI_lick_times = lick_start(ITI_lick_ind);
    padded_curr_ITI_lick_times = vertcat(0, curr_ITI_lick_times, 0);
    diff_lick = diff(padded_curr_ITI_lick_times);
    lick_bout_start_ind = find(diff_lick > 10);

    if ~isempty(lick_bout_start_ind)

         for j = 1:length(lick_bout_start_ind) 
             ITI_lick_time = curr_ITI_lick_times(lick_bout_start_ind(j));
             disp(lick_bout_start_ind(j))

            %% find frames for the trial

            fig_start = ITI_lick_time + t1; % start time for figure
            fig_end = ITI_lick_time + t2; % end time for figure

            ind_s = round(fig_start/dt);      %start index    
            ind_e = round(fig_end/dt);      %end index
            t_win = fig_start:dt:fig_end;


            %% find lick times
            [lick_ind,~] = find(lick_start > fig_start & lick_start < fig_end); %find licks during time win
            if ~isempty(lick_ind)
                lick_times = lick_start(lick_ind);
                lick_frames_rel = (lick_times - fig_start)*fr;  % lick frame RELATIVE to start of figure for current trial 
            end
            DF = df(:,ind_s:ind_e);     % raw dF/F for the image period only
            [r,n_frames] = size(DF);

            %% Get z-score
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
            space_h = round(r*.1);                               % number of rows
            space_holder = ones(space_h, n_frames)*5;
           
            if ~isempty(lick_ind)
                a = repelem(lick_frames_rel,2);             % x values
                y = [0 round(space_h*.9)];               
                b = repmat(y,1,length(lick_frames_rel));    % y values
                counter = 1;
            end    
           
            %% DF - z-scored 0 to 9
            figure(2)
            DF_z_sorted = vertcat(space_holder,DF_z_sorted);
            imagesc(DF_z_sorted);
            colormap(hot) 
            caxis([0 9]) 
            % step (3) - white patch 
            p1 = patch([0 n_frames+1 n_frames+1 0], [0 0 space_h+1 space_h+1] , 'white');
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
            colorbar
            title(strcat(mouse, {' '},day, {' ITI '}, num2str(trial)))
            % set(gca,'Visible','off')
            box off
            set(findall(gca, 'type', 'text'), 'visible', 'on')
            set(gcf, 'InvertHardCopy', 'off');
            %     save_path = strcat('C:\Users\clee162\Desktop\Figures\20-09-23\z-score 0 to 9\', celltype);
            %     s = strcat(save_path,filesep,celltype,'_',mouse,'_',day,'_',num2str(reward_num),'SORTED');
            %     if strcmp(save, 'Y') == 1
            %         saveas(gcf,s,'tif')
            %     end

            file_name_e = strcat(celltype,'_',mouse,'_ITI_',num2str(trial),'_heatmap.eps'); % .eps OR .tif
            file_name_t = strcat(celltype,'_',mouse,'_ITI_',num2str(trial),'_heatmap.tif'); % .eps OR .tif
            if strcmp(save, 'Y') == 1
                save_path = 'C:\Users\clee162\Desktop\Figures\20-11-13';
                s = strcat(save_path, filesep, file_name_e);
                print(gcf,'-depsc','-painters',s) % TO SAVE AS EPS
                s = strcat(save_path, filesep, file_name_t);
                saveas(gcf,s)   % TO SAVE AS TIFF
            end
        end
    end
end
