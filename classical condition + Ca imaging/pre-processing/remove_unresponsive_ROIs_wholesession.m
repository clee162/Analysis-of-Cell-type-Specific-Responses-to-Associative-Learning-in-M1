function real_rois = remove_unresponsive_ROIs_wholesession(df)
%last modified: 20-02-11

%find and exclude any ROIs that do not fire at least once in the session
%(irrespective of the behaviour)

%active ROIS defined by CHRIS HARVEY METHOD FROM DRISCOLL ET AL 2018

%shuffle matrix randomly 1000 times and find when
%dF/F is greater than the shuffled matrix 950 or more times. Function
%returns real_rois, a vector with the row # (ROI #) of active cells 

%% load data and set parameters
% df_filepath = 'D:\2P data\PV\CL147\Day 7\CL147_day7_df';
% ws_filepath = 'D:\2P data\PV\CL147\Day 7\CL147_day7_ws';
% load(df_filepath)
% load(ws_filepath)

[rois, frames] = size(df); 
threshold = 950;    %the number of times your df/f needs to be greater than the shuffled df/f
event_length = 5;   %this is the min number of frames considered to be an event

%% shuffle data
for i = 1:1000
    rand_shift = randi([1 frames]);
    df_shift = circshift(df, rand_shift,2); %randomly shift df/f
    
    temp_abv_shuffle = df > df_shift; %find where df/f is greater than the ith shift
    
    %sum the binary matrices to get the number of times each frame is above
    %the shuffled value
    if i == 1
        sum_abv_shuffle = temp_abv_shuffle; 
    else 
        sum_abv_shuffle = sum_abv_shuffle + temp_abv_shuffle; 
    end
end

%% Find runs equal or longer than event_length
for j = 1:rois
    abv_thresh = (sum_abv_shuffle(j,:) > threshold);  %where above threshold
    
    %aboveThreshold is a logical array, where 1 when above threshold, 0, below.
    %thus we want to calculate the difference between rising and falling edges
    abv_thresh = [false, abv_thresh, false];  %pad with 0's at ends
    edges = diff(abv_thresh);
    rising = find(edges==1);     %rising/falling edges
    falling = find(edges==-1);  
    spanwidth = falling - rising;  %width of span of 1's (above threshold)
    wideEnough = spanwidth >= event_length;   
    start_pos{j} = rising(wideEnough);    %start of each span
    end_pos{j} = falling(wideEnough)-1;   %end of each span
end

%% store row # of ROIs that were active at least once in real_rois
for roi = 1:numel(start_pos)
    if any(start_pos{roi})
                real_rois(roi) = roi;
    else 
        real_rois(roi) = 0;
    end 
    
end
