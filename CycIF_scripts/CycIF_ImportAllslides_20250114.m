%% CyCIF Import all slides (MCMICRO)
%  Jerry Lin 2025/01/14
%
%   Switch to the data directory first
%   Output: data tables & slideName
%

%% Initialization

disp('Please switch to the data directory first!!');
px1 = input('Please input pixel size (0.325 for Orion):');
sw1 = input('Remove zero (y/n)?','s');

%% Import data

dir1 = dir;
dir1 = struct2table(dir1);
dir1 = dir1.name;

filelist = dir1;
filelist = filelist(3:end);
slideName = cellfun(@(X) X(1:8),filelist,'UniformOutput',false);

tic;
for i = 1:length(slideName)
    filename = filelist{i};
    disp(strcat('Processing:',slideName{i},'.....',num2str(i)));
    data1 = CycIF_importMcMicro(filename,px1);
    if ismember(sw1,'y')
        data1 = CycIF_removezero(data1);
    end
    data1{:,2:end} = uint16(data1{:,2:end});
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

clear data1 dir1 filename i;
