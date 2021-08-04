%This runs through the folders for different animals and the subfolders for
%different days and uses DF traces to find unresponsive ROIs (these are
%ROIs that don't respond at least once during the entire session 
%This will save the file to the same folder as the DF 

clear 
clc

path = 'D:\2P data\VIP';
days = {'Day 1','Day 7'};
animals = {'CL201'};


for n = 1:length(animals)
    animal_path = strcat(path, filesep, animals{n});
    
    for day = 1:length(days)
        data_path = strcat(animal_path, filesep, days{day});
        matfiles = dir(fullfile(data_path,'*.mat'));
        df = [];
        for i = 1:length(matfiles) 
            s = (fullfile(data_path, matfiles(i).name));
            load(s) 
        end
        if ~isempty(df)
            real_rois = remove_unresponsive_ROIs_wholesession(df);

            file_name = strcat(animals{n},'_',days{day},'_responsive_ROI_indices');
            save_path = strcat(data_path, filesep, file_name);
            disp(save_path)
            save(save_path, 'real_rois')
            disp(strcat(file_name, ' saved'))
        else
            disp(strcat(animals{n}, days{day}, ' has no df file'))
        end
    end
    if n == length(animals) 
        disp('Finished')
    end
end
