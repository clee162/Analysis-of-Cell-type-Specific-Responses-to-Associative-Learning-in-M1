function [wt_start, lick_start, tone_start, tone_end, catch_t] = wavesurfer_analysis_2019_v4(ws_filepath)
%This function will output vectors with the water, lick and tone start times in seconds. 
%it removes catch trials from wt_start, tone_start and tone_end and stores
%them separately in catch_t. 

%wt_start - the times the water valve opens 
%lick_start - all the times the mouse breaks the beam break 
%tone_start - the time the tone starts, excluding catch trials 
%tone_end - the time the tone ends, excluding catch trials 
%catch_t - column 1: the trial #s that are catch trials . column 2: the
%tone_start time for catch trials

s = ws.loadDataFile(ws_filepath);
sample_rate = s.header.Acquisition.SampleRate;

%s.sweep_001.analogScans ia holds values for beam break in column 1, tone
%in column 2 and water in column 3 

if isfield(s,'sweep_0012')
     s.sweep_0001 = s.sweep_0012;
end 

raw_data.lick = s.sweep_0001.analogScans(:,1);
raw_data.tone = s.sweep_0001.analogScans(:,2);
raw_data.water = s.sweep_0001.analogScans(:,3); 

num_samples = length(s.sweep_0001.analogScans(:,1));

%make binary 
binary.lick = zeros(num_samples,1);
binary.tone = zeros(num_samples,1);
binary.water = zeros(num_samples,1);

for i = 1:num_samples
    if raw_data.lick(i) > 1
        binary.lick(i) = 1; 
    end
    
%     if raw_data.tone(i) > 7 
     if raw_data.tone(i) < 1 

        binary.tone(i) = 1; 
    end
    
    if raw_data.water(i) > 4.3 
        binary.water(i) = 1;
    end 
end

ind.lick = diff(binary.lick);
ind.tone = diff(binary.tone); 
ind.water = diff(binary.water); 

on_off_times.lick = (find(ind.lick))./sample_rate;
on_off_times.tone = (find(ind.tone))./sample_rate;
on_off_times.water = (find(ind.water))./sample_rate; 

wt_start = on_off_times.water(1:2:end);     %holds the start of water in seconds
lick_start = on_off_times.lick(1:2:end);    %holds the start of licks in seconds
tone_end_ind = find(diff([on_off_times.tone])>=60);
temp_end = on_off_times.tone((length(on_off_times.tone))); 
tone_end = on_off_times.tone(tone_end_ind);
tone_end = vertcat(tone_end,temp_end);
tone_start_ind = tone_end_ind + 1; %first value of jump
temp_start(1) = 1;
tone_start_ind = vertcat(temp_start, tone_start_ind);
tone_start = on_off_times.tone(tone_start_ind);  %holds the start of tone in seconds

%% in case of repeats in wt file 
d = diff(wt_start); 
m = find(d < 50);

for i = 1:length(m) 
    d = diff(wt_start); 
    m = find(d < 50);
    wt_start(m(1)+1) = [];
end

%% Remove catch trials and store them separately
wt_counter = 1;
tone_counter = 1;
catch_counter = 1;

for i = 1:length(tone_start) 
    d = (wt_start(wt_counter) - tone_start(tone_counter));
    if d < 2 || d > 3 % if d is less than 2 or greater than 3, this is not a normal trial
        catch_t(catch_counter,1) = i;
        catch_t(catch_counter,2) = tone_start(tone_counter);
        catch_counter = catch_counter + 1;
        tone_counter = tone_counter + 1;
    end
    
    if d > 2 && d < 3   % this is a regular trial
        wt_counter = wt_counter + 1;
        tone_counter = tone_counter + 1;
    end
end     
        
        

% if length(tone_start) == length(wt_start)
%     catch_t = [];
% end
% 
% if length(tone_start) > length(wt_start)
%     
%     num_catch = length(tone_start) - length(wt_start); 
%     counter = 1;
%     if num_catch > 0
%         catch_t = zeros(num_catch,2);
%         h = 1;
%         l = 1;
%         for i = 1:length(wt_start)
% 
%             if wt_start(l) - tone_start(l) > 10 %catch
%                 catch_t(h, 1) = i; 
%                 catch_t(h, 2) = tone_start(l); 
%                 tone_start(l) = [];
%                 h = h+1;
%                 l = l-1;
%             end
% 
%             while i == length(wt_start) && length(tone_start) > i
%                 catch_t(h, 1) = i + counter; 
%                 catch_t(h, 2) = tone_start(i+1); 
%                 tone_start(i) = [];
%                 counter = counter + 1; 
%                 h = h + 1; 
%             end
%             l = l+1;
%         end
%     end
   % tone_end(catch_t(:,1)) = [];
   
%    a = wt_start - tone_start;
%    
%    for i = 1:length(a) 
%        if a(i) > 2.7 || a(i) < 2 
%            wt_start(i) = tone_start(i) + 2.5;
%        end
%    end
           

% day = ws_filepath(end-4:end);
% animal = ws_filepath(end-10:end-6);
% new_filename = strcat(animal,'_',day,'_ws');
% save_path = strcat(ws_filepath,filesep,new_filename);
% save(save_path, 'wt_start', 'lick_start', 'tone_start', 'tone_end', 'catch_t')
% disp([save_path, ' saved']);
end