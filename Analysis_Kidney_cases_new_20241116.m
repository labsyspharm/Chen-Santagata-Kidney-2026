%% Analysis Kindey cases
%  Jerry Lin 2024/11/16


%% Import data & clean up

dir1 = dir;
dir1 = struct2table(dir1);
dir1 = dir1.name;

filelist = dir1;
filelist = filelist(3:end);
slideName = cellfun(@(X) X(1:8),filelist,'UniformOutput',false)

tic;
%for i = 1:length(slideName)
for i = 140:140
    filename = filelist{i};
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = CycIF_tif2table(filename,1.3);  %level 3 (Orion)
    data1 = data1(data1.ch1>exp(6.1),:);
    data1 = data1(1:4:end,:);
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Assign Omero ROIs

dir1 = dir;
dir1 = struct2table(dir1);
dir1 = dir1.name;

filelist = dir1;
ROIlist = filelist(3:end);

tic;
%for i = 1:length(slideName)
    i = 140;
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1 = CycIF_assignOmeroROIsingle2(data1,ROIlist{i},0.325,0,'GLOm',0);
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
%end

%% Reassign channel names
tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1.Properties.VariableNames(1:19) = label1;
    data1{:,:} = uint16(data1{:,:});

    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Clean up data
tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    %data1 = data1(data1.DNA>exp(6.1),:);
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Generate all Glo data (without z scores)

allGLOm2 = table;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1 = data1(:,labelT);
    data2 = data1(data1.GLOm>0,:);
    data2.slideName = repmat(slideName(i),size(data2,1),1);
    if isempty(allGLOm2)
        allGLOm2 = data2;
    else
        allGLOm2 = vertcat(allGLOm2,data2);
    end
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Scatter plots

figure,scatter(log(sumGLOm2.mean_Kappa),log(sumGLOm2.mean_Lambda),20,grp2idx(sumGLOm2.slideName),'fill');
colormap(lines);
xlabel('Kappa');
ylabel('Lambda');
title('Lambda/Kappa');
set(gcf,'color','w')
gname(sumGLOm2.slideName);
lsline

figure,scatter(log(sumGLOm2.mean_PODXL),log(sumGLOm2.mean_Synaptopodin),20,grp2idx(sumGLOm2.slideName),'fill');
colormap(lines);
xlabel('PODXL');
ylabel('Synaptopodin');
title('Synaptopodin/PODXL');
set(gcf,'color','w')
gname(sumGLOm2.slideName);
lsline

%% Heatmaps

figure, imagesc(test2(:,temp1))
colormap(redbluecmap);
colorbar;
caxis([-0.3 0.5]);
xlabel('Orion measurement');
ylabel('Diagnosis');
title('Glomeruli only');
set(gca,'xtick',1:18);
set(gca,'xticklabels',label2);
set(gca,'ytick',1:25);
set(gca,'yticklabels',labelX2(20:end));
set(gcf,'color','w')

figure, imagesc(test4(:,temp1))
colormap(redbluecmap);
colorbar;
caxis([-0.3 0.5]);
xlabel('Orion measurement');
ylabel('Diagnosis');
title('Whole sample');
set(gca,'xtick',1:18);
set(gca,'xticklabels',label2);
set(gca,'ytick',1:25);
set(gca,'yticklabels',labelX2(20:end));
set(gcf,'color','w')

%% All gates

allgates = table;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    %temp1 = double(data1{:,1:19});
    %temp1 = zscore(temp1);
    %data1{:,1:19} = temp1;
    data2 = data1(:,labelp);
    data2.slideName = repmat(slideName(i),size(data2,1),1);
    if isempty(allgates)
        allgates = data2;
    else
        allgates = vertcat(allgates,data2);
    end
    toc;
end

%% Generate allcorr
allcorr = zeros(length(labelp),length(labelp),length(slideName));

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    allcorr(:,:,i) = corr(data1{:,labelp});
end
toc;

%% plot heatmap

figure, imagesc(corr(allsample{allsample.GLOm>0,labelp}))
set(gca,'xtick',1:18);
set(gca,'ytick',1:18);
set(gca,'xticklabels',labelp);
set(gca,'yticklabels',labelp);
colormap(redbluecmap);
colorbar;

%% t-test by diagnosis

allpvalue = zeros(length(labelp),length(labelp));
diag1 = 'AcuteTubularInjury';
list1 = tableINFO{:,diag1};

tic;
for i = 1:length(labelp)
    for j = 1:length(labelp)
        allpvalue(i,j)=myttest2(squeeze(allcorr(i,j,:)),list1);
    end
end
toc;

figure,imagesc(allpvalue<0.01);
set(gca,'xtick',1:18);
set(gca,'ytick',1:18);
set(gca,'xticklabels',labelp);
set(gca,'yticklabels',labelp);
colormap(cool);
title(diag1,'FontSiz',16);

%% boxplot (allcorr/allpvalue);

%diag1 = 'DiffuseAndNodularDiabeticGlomerulosclerosis';
%list1 = tableINFO{:,diag1};

name1 = 'AQP1';
name2 = 'PODXL';

i = find(ismember(labelp,name1));
j = find(ismember(labelp,name2));

figure,myboxplot2d(squeeze(allcorr(i,j,:)),list1);
title(strcat(name1,'-',name2),'FontSize',16);

%% Generate sumAllsample;

tic;
for i = 1:20
    marker1 = strcat('topic',num2str(i));
    allsample{:,marker1}=false(size(allsample,1),1);
    allsample{allsample.topics==i,marker1} = true;
end
toc;

tic;
for i = 1:20
    marker1 = strcat('Kmean20_',num2str(i));
    allsample{:,marker1}=false(size(allsample,1),1);
    allsample{allsample.Kmean20==i,marker1} = true;
end
toc;

tic;
sumAllsample = varfun(@mean,allsample,'GroupingVariables','slideName');
sumAllsample = join(sumAllsample,tableINFO,'keys','slideName');
sumSampleGLOm = varfun(@mean,allsample(allsample.GLOm>0,:),'GroupingVariables','slideName');
sumSampleGLOm = join(sumSampleGLOm,tableINFO,'keys','slideName');
sumSampleOther = varfun(@mean,allsample(allsample.GLOm==0,:),'GroupingVariables','slideName');
sumSampleOther = join(sumSampleOther,tableINFO,'keys','slideName');
toc;

%% all markers boxplot by diagnosis

sum1 = sumSampleGLOm;
diag1 = 'AcuteInterstitialNephritis';

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:length(labelp)
    subplot(3,6,i);
    myboxplot2d(log(sum1{:,strcat('mean_',labelp{i})}),sum1{:,diag1});
    title(labelp{i});
end
sgtitle(strcat(diag1,'(Mean intensity)'),'FontSize',16);
set(gcf,'color','w')


%'AcuteInterstitialNephritis';
%'ChronicActiveInterstitialNephritis';

%% all topics boxplot by diagnosis

sum1 = sumSampleGLOm;
diag1 = 'IgArelatedGlomerularDisease';

figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:20
    subplot(3,7,i);
    myboxplot2d(sum1{:,strcat('mean_topic',num2str(i))}*100,sum1{:,diag1});
    ytickformat('percentage');
    title(strcat('topic',num2str(i)));
end
sgtitle(strcat(diag1,'(Mean intensity)'),'FontSize',16);
set(gcf,'color','w')


%'AcuteInterstitialNephritis';
%'ChronicActiveInterstitialNephritis';

%% Marker versus diagnosis

sum1 = sumAllsample;

arrayDiff = zeros(size(labelp,1),size(listDiag,1));
arrayPvalue = zeros(size(labelp,1),size(listDiag,1));
tic;
for i =1:length(labelp)
    for j=1:length(listDiag.Diag)
        marker1 = strcat('mean_',labelp{i});
        Diag1 = listDiag.Diag(j);
        [arrayPvalue(i,j),arrayDiff(i,j)] = myttest2(log(sum1{:,marker1}),sum1{:,Diag1});
    end
end
toc;

% Plot
arrayDiff = abs(arrayDiff);
arrayPvalue = -log(arrayPvalue);
arrayPvalue(isnan(arrayPvalue))=0;

xcorr = repmat((1:size(arrayDiff,1))',1,size(arrayDiff,2));
ycorr = repmat((1:size(arrayDiff,2)),size(arrayDiff,1),1);
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(xcorr(:),ycorr(:),150*arrayDiff(:)+20,arrayPvalue(:),'fill');
box off;

colormap(cool);
caxis([0 6.5]);
colorbar;

set(gcf,'color','w');
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp);
set(gca,'ytick',1:length(listDiag.label));
set(gca,'yticklabels',listDiag.label);


%%  Gated Marker versus diagnosis

sum1 = sumSampleGLOm;

arrayDiff = zeros(size(labelp,1),size(listDiag,1));
arrayPvalue = zeros(size(labelp,1),size(listDiag,1));
tic;
for i =1:length(labelp)
    for j=1:length(listDiag.Diag)
        marker1 = strcat('mean_',labelp{i});
        Diag1 = listDiag.Diag(j);
        [arrayPvalue(i,j),arrayDiff(i,j)] = myttest2(log(sum1{:,marker1}),sum1{:,Diag1});
    end
end
toc;

% Plot
arrayDiff = abs(arrayDiff);
arrayPvalue = -log(arrayPvalue);
arrayPvalue(isnan(arrayPvalue))=0;

xcorr = repmat((1:size(arrayDiff,1))',1,size(arrayDiff,2));
ycorr = repmat((1:size(arrayDiff,2)),size(arrayDiff,1),1);
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(xcorr(:),ycorr(:),200*arrayDiff(:)+20,arrayPvalue(:),'fill');
box off;

colormap(cool);
caxis([0 6.5]);
colorbar;

set(gcf,'color','w');
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp);
set(gca,'ytick',1:length(listDiag.label));
set(gca,'yticklabels',listDiag.label);


%%  Topics versus diagnosis

sum1 = sumAllsample;

arrayDiff = zeros(maxT,size(listDiag,1));
arrayPvalue = zeros(maxT,size(listDiag,1));

tic;
for i =1:maxT
    for j=1:length(listDiag.Diag)
        marker1 = strcat('mean_topic',num2str(i));
        disp(i);
        Diag1 = listDiag.Diag(j);
        [arrayPvalue(i,j),arrayDiff(i,j)] = myttest2(sum1{:,marker1},sum1{:,Diag1});
    end
end
toc;

% Plot
arrayDiff = abs(arrayDiff);
arrayPvalue = -log(arrayPvalue);
arrayPvalue(isnan(arrayPvalue))=0;

xcorr = repmat((1:size(arrayDiff,1))',1,size(arrayDiff,2));
ycorr = repmat((1:size(arrayDiff,2)),size(arrayDiff,1),1);
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(xcorr(:),ycorr(:),600*arrayDiff(:)+20,arrayPvalue(:),'fill');
box off;

colormap(cool);
caxis([0 6.5]);
colorbar;

set(gcf,'color','w');
set(gca,'xtick',1:maxT);
%set(gca,'xticklabels',labelp);
xlabel('Topics');
set(gca,'ytick',1:length(listDiag.label));
set(gca,'yticklabels',listDiag.label);

%% Plot LDA
figure
for topicIdx = 1:maxT
    subplot(3,ceil(maxT/3),topicIdx)
    temp1 = table;
    temp1.Word = labelp3;
    
    temp1.Count = ldaAll.TopicWordProbabilities(:,topicIdx);
    wordcloud(temp1,'Word','Count');
    title("Topic: " + topicIdx)
end


%%  Kmean versus diagnosis

sum1 = sumAllsample;

arrayDiff = zeros(maxT,size(listDiag,1));
arrayPvalue = zeros(maxT,size(listDiag,1));

tic;
for i =1:maxT
    for j=1:length(listDiag.Diag)
        marker1 = strcat('mean_Kmean20_',num2str(i));
        disp(i);
        Diag1 = listDiag.Diag(j);
        [arrayPvalue(i,j),arrayDiff(i,j)] = myttest2(sum1{:,marker1},sum1{:,Diag1});
    end
end
toc;

% Plot
arrayDiff = abs(arrayDiff);
arrayPvalue = -log(arrayPvalue);
arrayPvalue(isnan(arrayPvalue))=0;

xcorr = repmat((1:size(arrayDiff,1))',1,size(arrayDiff,2));
ycorr = repmat((1:size(arrayDiff,2)),size(arrayDiff,1),1);
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(xcorr(:),ycorr(:),600*arrayDiff(:)+20,arrayPvalue(:),'fill');
box off;

colormap(cool);
caxis([0 6.5]);
colorbar;

set(gcf,'color','w');
set(gca,'xtick',1:maxT);
%set(gca,'xticklabels',labelp);
xlabel('Cluster');
set(gca,'ytick',1:length(listDiag.label));
set(gca,'yticklabels',listDiag.label);

%% Change channel names;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1.Properties.VariableNames(51) = "TUBp";
    data1.Properties.VariableNames(52) = "VESp";
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Generate VESp

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1.VESp = data1.SMAp & ~data1.ColIIIp & data1.ColIVp;
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Generate TUBa

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1.TUBa = data1.TUBp>0 & data1.GLOa==0 & data1.VES2a==0;
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Generate Region

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1.Region = zeros(size(data1,1),1);
    data1.Region(data1.TUBa) = 1;
    data1.Region(data1.VES2a>0) = 2;
    data1.Region(data1.GLOm>0) = 3;
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Remove channels;

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1 = removevars(data1, "VESp");
    data1 = removevars(data1, "VESa");
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% display LDA

figure
for topicIdx = 1:maxT
    subplot(4,ceil(maxT/4),topicIdx)
    temp1 = table;
    temp1.Word = labelp3;
    
    temp1.Count = ldaAll.TopicWordProbabilities(:,topicIdx);
    wordcloud(temp1,'Word','Count');
    title("Topic: " + topicIdx)
end


%% Generate sumAllsample2

tic;
sumAllsample2 = varfun(@mean,allsample,'GroupingVariables','slideName');
sumAllsample2 = join(sumAllsample2,tableINFO,'keys','slideName');
sumSampleGLOm2 = varfun(@mean,allsample(allsample.Region==3,:),'GroupingVariables','slideName');
sumSampleGLOm2 = join(sumSampleGLOm2,tableINFO,'keys','slideName');
sumSampleTUB2 = varfun(@mean,allsample(allsample.Region==1,:),'GroupingVariables','slideName');
sumSampleTUB2 = join(sumSampleTUB2,tableINFO,'keys','slideName');
sumSampleVES2 = varfun(@mean,allsample(allsample.Region==2,:),'GroupingVariables','slideName');
sumSampleVES2 = join(sumSampleVES2,tableINFO,'keys','slideName');
sumSampleMISC2 = varfun(@mean,allsample(allsample.Region==0,:),'GroupingVariables','slideName');
sumSampleMISC2 = join(sumSampleMISC2,tableINFO,'keys','slideName');
toc;

%% Compare Variables between regions

list1 = {'sumSampleGLOm2','sumSampleVES2','sumSampleTUB2'};
label1 = {'Glo.','Ves.','Tub.'};

name1 = 'NMF9_2';
marker1 = strcat('mean_',name1);
array1 = [];

for i = 1:length(list1)
    sum1 = eval(strcat(list1{i}));
    temp1 = log(sum1{:,marker1}+5);
    if isempty(array1)
        array1 = temp1;
    else
        array1 = horzcat(array1,temp1);
    end
end

figure,myboxplot3d(array1);
set(gca,'xticklabels',label1);
title(name1);

%% Test all markers between regions

list1 = {'sumSampleGLOm2','sumSampleVES2','sumSampleTUB2'};
label1 = {'Glo.','Ves.','Tub.'};

figure('units','normalized','outerposition',[0 0 1 1]);

for j = 1:length(labelp)
    %name1 = 'AF';
    name1 = labelp{j};
    marker1 = strcat('mean_',name1);
    array1 = [];
    subplot(3,6,j);

    for i = 1:length(list1)
        sum1 = eval(strcat(list1{i}));
        temp1 = log(sum1{:,marker1}+5);
        if isempty(array1)
            array1 = temp1;
        else
            array1 = horzcat(array1,temp1);
        end
    end

    myboxplot3(array1);
    set(gca,'xticklabels',label1);
    title(name1);
end


%% Test all NMF between regions

list1 = {'sumSampleGLOm2','sumSampleVES2','sumSampleTUB2'};
label1 = {'Glo.','Ves.','Tub.'};

figure('units','normalized','outerposition',[0 0 1 1]);

for j = 1:9
    name1 = strcat('NMF9_',num2str(j));
    marker1 = strcat('mean_',name1);
    array1 = [];
    subplot(2,5,j);

    for i = 1:length(list1)
        sum1 = eval(list1{i});
        temp1 = log(sum1{:,marker1}+5);
        if isempty(array1)
            array1 = temp1;
        else
            array1 = horzcat(array1,temp1);
        end
    end

    myboxplot3(array1);
    set(gca,'xticklabels',label1);
    title(name1,'Interpreter','none');
end

set(gcf,'color','w');

%% tSNE scores

list1 = {'sumAllsample2','sumSampleGLOm2','sumSampleVES2','sumSampleTUB2'};
label1 = {'All','Glo.','Ves.','Tub.'};
t1 = zeros(223,2,length(list1));

tic;
for i = 1:length(list1)
    sum1 = eval(list1{i});
    temp1 = sum1{:,strcat('mean_',labelp)}; 
    t1(:,:,i) = tsne(temp1);
    toc;
end

%% tSNE plots

figure('units','normalized','outerposition',[0 0 1 1]);
label1 = {'Whole sample','Glomeruli','Vasculates','Tubules'};

diag1 = 25;
marker1 = listDiag.Diag(diag1);
name1 = listDiag.label(diag1);

for i = 1:length(list1)
    subplot(2,2,i);
    scatter(t1(:,1,i),t1(:,2,i),30,sum1{:,marker1},'fill');
    colormap(cool);
    title(label1{i});
    xlabel('tsne1');
    ylabel('tsne2');
end

sgtitle(name1);
set(gcf,'color','w');


%% tSNE plots (markers)

figure('units','normalized','outerposition',[0 0 1 1]);
label1 = {'Whole sample','Glomeruli','Vasculates','Tubules'};

p = 18;
marker1 = strcat('mean_',labelp{p});
name1 = labelp{p};

for i = 1:length(list1)
    subplot(2,2,i);
    scatter(t1(:,1,i),t1(:,2,i),25,log(sum1{:,marker1}),'fill');
    colormap(jet);
    colorbar;
    title(label1{i});
    xlabel('tsne1');
    ylabel('tsne2');
end

sgtitle(name1);
set(gcf,'color','w');


%%  NMF versus diagnosis

sum1 = sumAllsample2;

arrayDiff = zeros(9,size(listDiag,1));
arrayPvalue = zeros(9,size(listDiag,1));

tic;
for i =1:9
    for j=1:length(listDiag.Diag)
        marker1 = strcat('mean_NMF9_',num2str(i));
        disp(i);
        Diag1 = listDiag.Diag(j);
        [arrayPvalue(i,j),arrayDiff(i,j)] = myttest2(sum1{:,marker1},sum1{:,Diag1});
    end
end
toc;

% Plot
arrayDiff = abs(arrayDiff);
arrayPvalue = -log(arrayPvalue);
arrayPvalue(isnan(arrayPvalue))=0;

xcorr = repmat((1:size(arrayDiff,1))',1,size(arrayDiff,2));
ycorr = repmat((1:size(arrayDiff,2)),size(arrayDiff,1),1);


figure('units','normalized','outerposition',[0 0 1 1]);
scatter(xcorr(:),ycorr(:),0.5*arrayDiff(:)+5,arrayPvalue(:),'fill');
box off;

colormap(cool);
caxis([0 6.5]);
colorbar;

set(gcf,'color','w');
set(gca,'xtick',1:9);
%set(gca,'xticklabels',labelp);
xlabel('NMF vectors');
xlim([0 10]);

set(gca,'ytick',1:length(listDiag.label));
set(gca,'yticklabels',listDiag.label);


%% all topics boxplot by diagnosis

sum1 = sumAllsample2;
%diag1 = 'DiffusePodocytopathy';
i = 5;
diag1 = listDiag.Diag{i};
name1 = listDiag.label{i};

figure('units','normalized','outerposition',[0.5 0 0.5 1]);

for i = 1:9
    subplot(3,3,i);
    myboxplot2(log(sum1{:,strcat('mean_NMF9_',num2str(i))}),sum1{:,diag1});
    title(strcat('NMF9_',num2str(i)),'Interpreter','none','FontSize',9);
end
sgtitle(name1,'FontSize',14);
set(gcf,'color','w')

%% all markers boxplot by diagnosis

sum1 = sumAllsample2;
%diag1 = 'DiffusePodocytopathy';
i = 5;
diag1 = listDiag.Diag{i};
name1 = listDiag.label{i};

figure('units','normalized','outerposition',[0.5 0 0.5 1]);

for i = 1:18
    subplot(4,5,i);
    myboxplot2(log(sum1{:,strcat('mean_',labelp{i})}),sum1{:,diag1});
    title(labelp{i},'Interpreter','none','FontSize',9);
end
sgtitle(name1,'FontSize',14);
set(gcf,'color','w')

%% PLS regression

figure;
stat1 = 'mean';
target1 = 'Nephrectomy';

%list1 = vertcat(strcat('mean_',labelp),strcat('std_',labelp));
%X = sumGLO_all{:,list1};

X = sumGLO_all{:,strcat(stat1,'_',labelp)};
y = sumGLO_all{:,target1};
[n,p] = size(X);
%p = 10;

[Xloadings,Yloadings,Xscores,Yscores,betaPLS,PLSPctVar] = plsregress(X,y,p);
subplot(1,2,1)
plot(1:p,cumsum(100*PLSPctVar(2,:)),'-bo');
ylims = ylim;
ylim([0,ylims(2)]);

ylabel ('Percent Variance explained');
xlabel ('Number of PLS components');
title(strcat(target1,'(',stat1,')'),'FontSize',20,'interpreter','none');

subplot(1,2,2)
imagesc(Xloadings);
colormap(redbluecmap);
set(gca,'ytick',1:p);
set(gca,'yticklabels',labelp);
title('PLS loadings');
xlabel('PLS components');


%% mean versus std

figure,scatter(log(sumGLO_all.mean_C1q),log(sumGLO_all.std_C1q),25,sumGLO_all.C1q_score,'fill');
colormap(jet);
xlabel('mean','FontSize',16);
ylabel('std','FontSize',16);
title('C1q','FontSize',18);

%% boplot (mean versus std)

name1 = 'Albumin';

figure,boxplot(log(sumGLO_all{:,strcat('mean_',name1)}),sumGLO_all{:,strcat(name1,'_score')});
corr(log(sumGLO_all{:,strcat('mean_',name1)}),sumGLO_all{:,strcat(name1,'_score')})
ylabel('log(intensity)');
xlabel('Clinical scores');
title(strcat(name1,' (mean)'),'FontSize',20);

figure,boxplot(log(sumGLO_all{:,strcat('std_',name1)}),sumGLO_all{:,strcat(name1,'_score')});
corr(log(sumGLO_all{:,strcat('std_',name1)}),sumGLO_all{:,strcat(name1,'_score')})
ylabel('log(intensity)');
xlabel('Clinical scores');
title(strcat(name1,' (std)'),'FontSize',20);

%% Import Yu-An's mask

tic;
for i = 1:length(slideName)

    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    filename = strcat(slideName{i},'-tubule_type_mask.ome.tif');
    temp1 = CycIF_tif2table(filename,1.3);
    temp1 = temp1(temp1.ch1>0,:);
    data1.Mask = zeros(size(data1,1),1);
    [idx,d]=knnsearch([temp1.Xt,temp1.Yt],[data1.Xt,data1.Yt]);
    list1 = temp1.ch1(idx);
    data1.Mask(d<10)= list1(d<10);
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%--- Add VESm column----

tic;
for i = 1:length(slideName)
    disp(strcat('Processing:',slideName{i},'....',num2str(i)));
    data1 = eval(strcat('data',slideName{i}));
    data1.VESm = data1.Mask==2;
    eval(strcat('data',slideName{i},'=data1;'));
    toc;
end

%% Compartment scatter plot
name1 = 'IgG';
marker1 = strcat('mean_',name1);

figure,scatter(log(sumAll_GLO{:,marker1}),log(sumAll_TUB{:,marker1}),40,sumAll{:,79},'fill');
colormap(jet);
hcb = colorbar;
xlabel('Glomeruli');
ylabel('Tubules');
title(hcb,'Vessels');
title(name1,'FontSize',20);

%% Compartment scatter plot (diagnosis)

figure('units','normalized','outerposition',[0 0 1 1]);
name1 = 'AF';
marker1 = strcat('mean_',name1);

for i = 1:25
    diag1 = listDiag.Diag(i);
    subplot(5,5,i)
    gscatter(log(sumAll_TUB{:,marker1}),log(sumAll_VES{:,marker1}),sumAll{:,diag1},'br','..',10);
    %colormap(jet);
    %hcb = colorbar;
    %xlabel(strcat(name1,'(Glomeruli)'));
    %ylabel(strcat(name1,'(Tubules)'));
    
    title(listDiag.label(i),'FontSize',8);
    legend off;
end
legend('no','yes');
set(gcf,'color','w');

%% Compartment scatter plot (diagnosis)

name1 = 'IgA';
marker1 = strcat('mean_',name1);
i = 25;
diag1 = listDiag.Diag(i);

figure;
gscatter(log(sumAll_GLO{:,marker1}),log(sumAll_TUB{:,marker1}),sumAll{:,diag1},'br','..',20);
%colormap(jet);
%hcb = colorbar;
xlabel(strcat(name1,'(Glomeruli)'));
ylabel(strcat(name1,'(Tubules)'));
title(listDiag.label(i),'FontSize',16);
legend('No','Yes');

%% Bi-scatter plot (diagnosis)

figure('units','normalized','outerposition',[0 0 1 1]);
name1 = 'ColIII';
name2 = 'ColIV';
marker1 = strcat('mean_',name1);
marker2 = strcat('mean_',name2);

for i = 1:25
    diag1 = listDiag.Diag(i);
    subplot(5,5,i)
    gscatter(log(sumAll_GLO{:,marker1}),log(sumAll_GLO{:,marker2}),sumAll{:,diag1},'br','..',10);
    %colormap(jet);
    %hcb = colorbar;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    xlabel(name1);
    ylabel(name2);
    
    title(listDiag.label(i),'FontSize',8);
    legend off;
end
legend('no','yes');
set(gcf,'color','w');

%% calculate all cluster score

arrayClusterScore = zeros(length(labelp),length(labelp),length(slideName));

tic;
for s = 153:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    for i = 1:length(labelp)
        marker1 = labelp{i};
        for j = 1:length(labelp)
            marker2 = labelp{j};
           
            gate1 = strcat(marker1,'p');
            gate2 = strcat(marker2,'p');
            
            data2 = data1(data1{:,gate1},:);
            data3 = data1(data1{:,gate2},:);
            data4 = datasample(data1,size(data3,1));
            if size(data2,1)>0 & size(data3,1)>0
                [~,d1] = knnsearch([data2.Xt,data2.Yt],[data3.Xt,data3.Yt],'K',2);
                [~,d2] = knnsearch([data2.Xt,data2.Yt],[data4.Xt,data4.Yt],'K',2);
                arrayClusterScore(i,j,s) = mean(d2(:,2))/mean(d1(:,2));
            end
        end
    end
    toc;
end

%% check cluster score per sample

i = 120;

figure,imagesc(arrayClusterScore(:,:,i));
colormap(jet);
colorbar;

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
caxis([0 5]);
title(slideName{i},'FontSize',18);
set(gcf,'color','w');

%% Test all pvalues

arrayClusterPvalues = zeros(length(labelp),length(labelp));
diag1 = 'AcuteTubularInjury';

tic;
for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayClusterPvalues(i,j)=myttest2(squeeze(arrayClusterScore(i,j,:)),sumAll{:,diag1});
    end
end
toc;

figure,imagesc(arrayClusterPvalues<0.001);

set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);

temp1 = arrayClusterPvalues<0.001;
figure,circularGraph(temp1,'Label',labelp);
%% check 

name1 = 'CD8a';
name2 = 'SMA';

i = find(ismember(labelp,name1));
j = find(ismember(labelp,name2));

figure,myboxplot2d(squeeze(arrayClusterScore(i,j,:)),sumAll.Response);
title(strcat(name1,'-',name2),'FontSize',18);

%% Check auto-correlation and cluster side-by-side

i = 150;
figure('units','normalized','outerposition',[0 0 1 1]);

data1 = eval(strcat('data',slideName{i}));

subplot(1,2,1)
imagesc(corr(data1{:,labelp}));
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp);
title('Auto-correlation','FontSize',14);
colormap(jet);
colorbar('southoutside');
caxis([-0.25,0.75]);

subplot(1,2,2)
imagesc(arrayClusterScore(:,:,i));
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp);
title('Clustering Score','FontSize',14);
colormap(jet);
colorbar('southoutside');
caxis([0,5]);
set(gcf,'color','w');

sgtitle(slideName{i},'FontSize',18)


%% Check p-value side-by-side

diag1 = 'DiffusePodocytopathy';
cutoff1 = 0.01;
figure('units','normalized','outerposition',[0 0 1 1]);
tic;
% ----- Correlation -----
subplot(1,2,1);

arrayCorrPvalues = zeros(length(labelp),length(labelp));
for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayCorrPvalues(i,j)=myttest2(squeeze(allcorr(i,j,:)),sumAll{:,diag1});
    end
end
toc;

imagesc(arrayCorrPvalues<cutoff1);
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
title('Auto-correlation','FontSize',14)

% ---- Cluster ----
subplot(1,2,2)

arrayClusterPvalues = zeros(length(labelp),length(labelp));
for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayClusterPvalues(i,j)=myttest2(squeeze(arrayClusterScore(i,j,:)),sumAll{:,diag1});
    end
end
toc;

imagesc(arrayClusterPvalues<cutoff1);
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
title('Close-Proximity','FontSize',14)
set(gcf,'color','w');
sgtitle(diag1,'FontSize',18)

temp1 = arrayClusterPvalues<cutoff1 & ~(arrayCorrPvalues<cutoff1);
figure,circularGraph(temp1,'Label',labelp);
sgtitle(diag1,'FontSize',18);

%% Check p-value side-by-side (Nephrectomy as control)

ref1 = 'Nephrectomy';
%ref1 = 'AcuteInterstitialNephritis';
diag1 = 'DiffusePodocytopathy';
cutoff1 = 0.001;
figure('units','normalized','outerposition',[0 0 1 1]);

flag1 = sumAll{:,ref1}==1;
flag2 = sumAll{:,diag1}==1;
list1 = repmat({ref1},sum(flag1),1);
list2 = repmat({diag1},sum(flag2),1);
listAll = vertcat(list1,list2);

% ----- Correlation -----
subplot(1,2,1);
arrayCorrPvalues = zeros(length(labelp),length(labelp));

arrayTemp1 = allcorr(:,:,flag1);
arrayTemp2 = allcorr(:,:,flag2);
arrayTempAll = cat(3,arrayTemp1,arrayTemp2);

for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayCorrPvalues(i,j)=myttest2(squeeze(arrayTempAll(i,j,:)),listAll);
    end
end
toc;

imagesc(arrayCorrPvalues<cutoff1);
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
title('Auto-correlation','FontSize',14)

% ---- Cluster ----
subplot(1,2,2)

arrayTemp1 = arrayClusterScore(:,:,flag1);
arrayTemp2 = arrayClusterScore(:,:,flag2);
arrayTempAll = cat(3,arrayTemp1,arrayTemp2);

arrayClusterPvalues = zeros(length(labelp),length(labelp));
for i = 1:length(labelp)
    for j = 1:length(labelp)
        arrayClusterPvalues(i,j)=myttest2(squeeze(arrayTempAll(i,j,:)),listAll);
    end
end
toc;

imagesc(arrayClusterPvalues<cutoff1);
set(gca,'xtick',1:length(labelp));
set(gca,'xticklabels',labelp3);
set(gca,'ytick',1:length(labelp));
set(gca,'yticklabels',labelp3);
title('Close-Proximity','FontSize',14)
set(gcf,'color','w');
sgtitle(strcat(diag1,{' vs '},ref1),'FontSize',18)

tempInt = arrayClusterPvalues<cutoff1 & ~(arrayCorrPvalues<cutoff1);
figure,circularGraph(tempInt,'Label',labelp);
sgtitle(diag1,'FontSize',18);

%% Calculate distance to GLOm

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    data2 = data1(data1.GLOm>0,:);
    [~,data1.dGLO]=knnsearch([data2.Xt,data2.Yt],[data1.Xt,data1.Yt],'k',1);
    eval(strcat('data',slideName{s},'=data1;'));
    toc;
end

%% check 

name2 = 'KIM1';
name1 = 'Synaptopodin';

i = find(ismember(labelp,name1));
j = find(ismember(labelp,name2));

figure,
subplot(1,2,1);
myboxplot2d(squeeze(allcorr(i,j,:)),sumAll{:,diag1});
xlabel(diag1);
%set(gca,'xticklabels',{'-','+'});
ylabel('auto-correlation');

subplot(1,2,2);
myboxplot2d(squeeze(arrayClusterScore(i,j,:)),sumAll{:,diag1});
xlabel(diag1);
%set(gca,'xticklabels',{'-','+'});
ylabel('Clustering scores');

sgtitle(strcat(name1,'-',name2),'FontSize',18); 

%% check all KIM1 interaction

diag1 = 'DiffusePodocytopathy';

figure('units','normalized','outerposition',[0 0 1 1]);

name2 = 'KIM1';
for i = 1:length(labelp)
    subplot(3,6,i)
    name1 = labelp{i};
    %i = find(ismember(labelp,name1));
    j = find(ismember(labelp,name2));

    myboxplot2d(squeeze(allcorr(i,j,:)),sumAll{:,diag1});
    xlabel(diag1);
    %set(gca,'xticklabels',{'-','+'});
    ylabel('Clustering scores');
    title(strcat(name1,'-',name2),'FontSize',11); 
end

sgtitle(diag1,'FontSize',16); 

%% All generate sumKIM1c 

sumKIM1c = table;

tic;
for s = 1:length(slideName)
    name1 = slideName{s};
    disp(strcat('Processing:',name1,'.....',num2str(s)))
    data1 = eval(strcat('data',name1));

    sum1 = varfun(@mean,data1,'GroupingVariables','KIM1c');
    sum1.slideName = repmat(slideName(s),size(sum1,1),1);
    if isempty(sumKIM1c)
        sumKIM1c = sum1;
    else
        sumKIM1c = vertcat(sumKIM1c,sum1);
    end
    toc;
end

%% Check KIM1 size & cluster

diag1 = 'DiffusePodocytopathy';

figure,
subplot(1,2,1);
myboxplot2d(sumKIM1c2.GroupCount./sumAll.GroupCount *50000,sumAll{:,diag1});
xlabel(diag1);
%set(gca,'xticklabels',{'-','+'});
title('Cluster counts');

subplot(1,2,2);
myboxplot2d(log(sumKIM1c2.mean_GroupCount),sumAll{:,diag1});
xlabel(diag1);
%set(gca,'xticklabels',{'-','+'});
title('Cluster size');

%sgtitle(diag1,'FontSize',18); 

%% Analysis Orion & TITAN feature tables

%load('matlab-Kindey_cases_summary-20250520.mat')
% Col 1-768 Titan features
% Col 769-848 Orion features
% Col 862-886 Diagnosis

sum1 = tableCombine;

labelX = sum1.Properties.VariableNames;
labelX = labelX';
labelTitan = labelX(1:768);  % 768
labelOrion = labelX(769:848);  % 80
labelDiag = labelX(862:886);   % 25

countDiag = sum(sum1{:,labelDiag});
countDiag = countDiag';

arrayPvalueTitan = zeros(length(labelDiag),length(labelTitan));

tic;
for i =1:length(labelDiag)
    for j= 1:length(labelTitan)
        arrayPvalueTitan(i,j) = myttest2(sum1{:,labelTitan{j}},sum1{:,labelDiag{i}});
    end
end
toc;

arrayPvalueOrion = zeros(length(labelDiag),length(labelOrion));

tic;
for i =1:length(labelDiag)
    for j= 1:length(labelOrion)
        arrayPvalueOrion(i,j) = myttest2(sum1{:,labelOrion{j}},sum1{:,labelDiag{i}});
    end
end
toc;



%% Output plot (Titan Pvalue)

temp1= arrayfun(@(X) {strcat('(',num2str(X),')')},countDiag);
temp2 = strcat(labelDiag,temp1);

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,12,1:11);
imagesc(arrayPvalueTitan<0.01);
set(gca,'ytick',1:length(labelDiag));
set(gca,'yticklabels',temp2);
xlabel('Titan feature');
set(gca,'xtick',[]);

subplot(1,12,12)
barh(flip(sum(arrayPvalueTitan<0.01,2)));
set(gca,'ytick',[]);
title ('Counts(p<0.01)');
ylim([0.5 length(labelDiag)+0.5]);

sgtitle('Titan features versus Diagnosis','FontSize',18);
set(gcf,'color','w');

%% Output plot (Orion Pvalue)

cutoff = 0.01;

temp1= arrayfun(@(X) {strcat('(',num2str(X),')')},countDiag);
temp2 = strcat(labelDiag,temp1);

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,12,1:11);
imagesc(arrayPvalueOrion<cutoff);
set(gca,'ytick',1:length(labelDiag));
set(gca,'yticklabels',temp2);
%xlabel('Orion feature');
set(gca,'xtick',[]);
ylims = ylim;
line([20.5 20.5],ylims,'Color','k','LineWidth',1.5);
line([40.5 40.5],ylims,'Color','k','LineWidth',1.5);
line([60.5 60.5],ylims,'Color','k','LineWidth',1.5);
text(8,26,'GLO features','FontSize',12);
text(28,26,'INT features','FontSize',12);
text(48,26,'TUB features','FontSize',12);
text(68,26,'VES features','FontSize',12);

subplot(1,12,12)
barh(flip(sum(arrayPvalueOrion<cutoff,2)));
set(gca,'ytick',[]);
title ('Counts(p<0.01)');
ylim([0.5 length(labelDiag)+0.5]);

sgtitle('Orion features versus Diagnosis','FontSize',18);
set(gcf,'color','w');




