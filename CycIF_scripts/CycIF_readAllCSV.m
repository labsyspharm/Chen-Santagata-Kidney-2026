%% Import All csv files (MCMICRO)
%  Jerry Lin  2202/11/05
%
%  chdir to the csv file location
%


%% --- Generate filelist and slideName ---
test1 = dir;
test1 = struct2table(test1);
test1 = test1.name;
filelist = test1(3:end);
slideName = strrep(filelist,'.csv','');

%% --- Read all csv ----
tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data1 = CycIF_importMcMicro(filelist{i},0.65);
    data1{:,2:end} = uint16(data1{:,2:end});
    data1 = CycIF_removezero(data1);
%     names = data1.Properties.VariableNames;
%     idx1 = find(ismember(names,'AREA'))-1;
%     if(isempty(idx1))
%         idx1 = find(ismember(names,'Area'))-1;
%     end
%     idx1 = idx1 -6;
%     data1 = data1(data1{:,idx1}>exp(7.5),:);
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

