%%  CycIF_Generalsegmentation_20240804
%
%   Segmentation positive markers for all slides
%   Jerry Lin 2024/08/04

%% Initialization

gate1 = input('Please input the segmentated gate:','s');  
output1 = input('Please input the output name:','s');

minD = input('Please input minimal distance (default=50):');   %100;
minSize = input('Please input minimal cell count (default=5):'); %100;
flag1 = input('Output maps (true/false):');

%% Processing all slides

tic;
for i =1:length(slideName)
    disp(strcat('Processing:',slideName{i}));
    data0 = eval(strcat('data',slideName{i}));
    
    if ~ismember(data0.Properties.VariableNames,'CellID')
        data0.CellID = (1:size(data0,1))';
    end
    
    data1 = data0(data0{:,gate1}>0,:);
    
    if size(data1,1) > 0 
        pc1 = pointCloud(double([data1.Xt,data1.Yt,ones(size(data1,1),1)]));
        label1=pcsegdist(pc1,minD);
        data1.label1 = label1;
    
        table1 = tabulate(data1.label1);
        table1 = table1(table1(:,2)> minSize,:);
        data1 = data1(ismember(data1.label1,table1(:,1)),:);
        table1(:,4) = (1:size(table1,1))';
        data1.label2 = zeros(size(data1,1),1);
        for j = 1:size(table1,1)
            data1.label2(data1.label1==table1(j,1)) = table1(j,4);
        end
    
        data0{:,output1} = zeros(size(data0,1),1);
        data0{ismember(data0.CellID,data1.CellID),output1} = data1.label2;
        eval(strcat('data',slideName{i},'=data0;'));
        tabulate(data0{:,output1});
        toc;
    else
        data0{:,output1} = zeros(size(data0,1),1);
        eval(strcat('data',slideName{i},'=data0;'));
        tabulate(data0{:,output1});
        toc;
    end
    %clear data0 data1;
end

%% Output Segmeanted maps

if flag1 
    tic;
    for i =1:length(slideName)
        disp(strcat('Processing:',slideName{i}));
        data1 = eval(strcat('data',slideName{i}));
    
        figure('units','normalized','outerposition',[0 0 1 1]);
    
        CycIF_tumorview(data1,output1,9,100000);
        daspect([1 1 1]);
        set(gcf,'color','w');
        title(slideName{i});
        filename = strcat(slideName{i},'_',output1,'.png');
        saveas(gcf,filename);
        close;
        toc;
    end
end  %if
