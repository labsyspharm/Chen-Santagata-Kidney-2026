%% Apply ROI to all dataset
% 2025/11/17  Jerry Lin
%

disp('Please first switch to the folder contains all the ROI files!!');
flag1 = input('Keep the ROI name column? (true/false):');
scale1 = input('Scale facotr (0.325 for Orion):');
dir1 = dir;
dir1 = struct2table(dir1);
listROI = dir1.name;
listROI = listROI(3:end);

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data1 = eval(strcat('data',slideName{i}));
    data1.ROI = zeros(size(data1,1),1);
    %data1.ROIname = repmat({'none'},size(data1,1),1);
    filename = listROI{i};
    if exist(filename,'file')
        data1 = CycIF_assignOmeroROIsingle2(data1,filename,scale1,flag1,'ROI',true);
    end
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

disp('Completed!');
