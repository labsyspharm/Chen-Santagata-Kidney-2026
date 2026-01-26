%% CyCIF LDA count table
%  Jerry Lin 2023/06/24
%  
%  Require variables:  slideName, labelp2
%  Output:  allcounts

%% Counts tables for all slides (only markers)

tic;
maxL = input('Please input grid distance(default = 200um):');
allcounts = [];

for s=1:length(slideName)
    disp(strcat('Processing:',slideName{s}));
    data1 = eval(strcat('data',slideName{s}));

    data1.col = round(data1.Xt ./ maxL)+1;
    data1.row = round(data1.Yt ./ maxL)+1;
    data1.frame = round(data1.col+ max(data1.col)*(data1.row-1));

    count1 = varfun(@sum,data1,'GroupingVariables','frame','Inputvariables',labelp2);
    if isempty(allcounts)
        allcounts = count1{:,3:end};
    else
        allcounts = vertcat(allcounts,count1{:,3:end});
    end
    eval(strcat('data',slideName{s},'=data1;'));
    toc;
end   