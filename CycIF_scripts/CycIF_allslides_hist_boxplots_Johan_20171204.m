%% Plot all samples with histograms & boxplot
%  Jerry Lin 2017/12/04

slideName = {'Parental_T2',...
'Parental_T3',...
'Parental_T4',...
'SCP01B_T1',...
'SCP01B_T2',...
'SCP01B_T3',...
'SCP3B_T1',...
'SCP3B_T2',...
'SCP3B_T3',...
'SCP17B_T1',...
'SCP17B_T2',...
'SCP17B_T3',...
'SCP29_T1',...
'SCP29_T2',...
'SCP29_T3',...
'SCP32_T1',...
'SCP32_T2',...
'SCP32_T3'};

marker1 = 'Ki67';
population = 'EGFR';
FDR = 0.01;

eval(strcat('temp1=',population,slideName{1},';'));
temp2 = NaN(length(temp1.X),length(slideName));
temp3 = NaN(length(slideName),1);

%% Global gating
data1 = log(allsample{:,marker1});

figure,[gate, pluscells]=findcutoff(data1,2,FDR);
[~,~,~,lowb,highb]=findgate3(data1,0,0.05,0);
xlim([lowb-0.5,gate+3]);
title (strcat('Allsample gating:',marker1));

%% Indivual histogram
figure;

for slide = 1:length(slideName)
    myName = slideName{slide};
    eval(strcat('temp1=',population,myName,';'));
    data1 = log(temp1{:,marker1});
    subplot(3,6,slide);

    %[pluscells, gate, peak,lowb,highb] = findgate3(data1,1,0.05,0);
    [pluscells, gate, ~,lowb,highb]=findgate3(data1,1,0.05,gate);
    %[gate, pluscells]=findcutoff(data1,2,0.01);
    legend off;
    xlim([lowb-0.5,gate+3]);
    temp3(slide) = mean(pluscells);
    
    title(slideName(slide));
    ylabel(strcat(marker1,{' Intensity'}));
    temp2(1:length(data1),slide) = data1;
end
    
    
%% boxplot for all samples;

figure,boxplot(temp2,'Notch','on');
set(gca,'xtick',1:length(slideName));
set(gca,'xticklabels',slideName);
set(gca,'XTickLabelRotation',90);
ylabel(strcat(marker1,{' Intensity'}));

%% plot positive cells

temp3 = reshape(temp3,3,6);
figure,bar(temp3');
names = {'Parental','SCP01B','SCP3B','SCP17B','SCP29','SCP32'};
set(gca,'xticklabels',names);
set(gca,'XTickLabelRotation',45);
ylabel(strcat(marker1,{' Postive'}));