% plot total number of responsive ROIs by CH method irrespective of
% behaviour 

% This plots the total percent of active cells within a session for PN, PV,
% VIP and SOM (Fig 2D) 

clear 
close all
clc
%% Enter info for group you want to analyze
PN_mice_IDs ={'CL174', 'CL175', 'CL176', 'CL181', 'CL182', 'CL184'};
% % 
PV_mice_IDs = {'CL147','CL172', 'CL195', 'CL196', 'CL198','CL158'};
% % 
VIP_mice_IDs = {'CL136', 'CL144', 'CL146', 'CL188'};

SOM_mice_IDs = {'CL170', 'CL173', 'CL186', 'CL187', 'CL191', 'CL192', 'CL193'};


day = 'Day 7';
outer_folder = 'D:\2P data';
save = 'n'; %save = 'Y' if you want to save
save_path = 'C:\Users\clee162\Desktop\Figures\21-06-11';

%% Now start analysis
celltypes = {'PN','PV','VIP','SOM'};
mice_IDs{1} = PN_mice_IDs;
mice_IDs{2} = PV_mice_IDs;
mice_IDs{3} = VIP_mice_IDs;
mice_IDs{4} = SOM_mice_IDs;

for celltype = 1:length(celltypes)
    curr_celltype_string = celltypes{celltype}; 
    curr_celltype_mice_IDs = mice_IDs{celltype};
    n = length(curr_celltype_mice_IDs);
    
    for mouse = 1:n
        current_mouse = curr_celltype_mice_IDs{mouse};
        path = strcat(outer_folder, filesep, curr_celltype_string, filesep, current_mouse, filesep, day);
        disp(path)
        matfiles =  dir(fullfile(path, '*.mat'));
        nfiles = length(matfiles);
            for j = 1:length(matfiles)
                s = (fullfile(path, matfiles(j).name));
                load(s)
            end
            
            [total_rois,~] = size(df);
            real_rois(real_rois == 0) = [];
            number_active_rois = length(real_rois);
            percent_active_rois{celltype}(mouse) = (number_active_rois/total_rois)*100;
    end
    
    overall_avg(celltype) = mean(percent_active_rois{celltype});
    sem(celltype) = std(percent_active_rois{celltype})/sqrt(length(percent_active_rois{celltype}));
    
end

data_for_anova = cell2mat(percent_active_rois);
groupings = {'PN','PN','PN','PN','PN','PN','PV','PV','PV','PV','PV','PV','VIP','VIP','VIP','VIP','SOM','SOM','SOM','SOM','SOM','SOM','SOM',};
[p,~,~] = anova1(data_for_anova,groupings,[],'off');

a = [0 0.12 -0.1 -0.02 0.04 -0.09 0.12 -0.05 0.03];
x = 1:1:length(celltypes);

for j = 1:length(celltypes)
    disp(strcat(day, {' '}, celltypes(j), {': '}, num2str(round(overall_avg(j),2)), {' ± '}, num2str(round(sem(j),2))));
end


%% figure
close all
figure(1)
bar(x,overall_avg,'FaceColor',[1 1 1])
hold on 
er = errorbar(x,overall_avg,sem, sem);
er.Color = [0 0 0];
er.LineStyle = 'none';  
hold on
sz = 55;
c(1,:) = [0 0 0];
c(2,:) = [0.9 0 0];
c(3,:) = [0 0.5 1];
c(4,:) = [0.3 0.7 0.25];
for celltype = 1:length(celltypes)
    for point = 1:length(mice_IDs{celltype})
        scatter_x = a(point)+celltype;
        scatter(scatter_x,percent_active_rois{celltype}(point),sz,c(celltype,:),'filled',...
            'MarkerFaceAlpha', 0.8,...
            'MarkerEdgeColor', 'None')
        hold on
    end
end
xlim([0 length(celltypes)+1])
ylim([0 120])
ylabel('Percent of ROIs')
title('Overall active ROIs - whole session (CH method)')
xticks(1:length(celltypes))
% yticks(0:10:100)
xticklabels(celltypes)
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1), NE(2), ['ANOVA: p =',num2str(p)], 'VerticalAlignment','top', 'HorizontalAlignment','right');
if strcmp(save,'Y') == 1
    file_name = 'percent_active_whole_session';
    s = strcat(save_path, filesep, file_name);
    print(gcf,'-depsc','-painters',s)
end
