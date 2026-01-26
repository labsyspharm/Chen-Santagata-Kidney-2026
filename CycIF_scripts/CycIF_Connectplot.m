function [h1,h2]=CycIF_Connectplot(data1,labelp,k,sw1)
%% For viewing a specific channel in 'real space', require CycIF data table
%  Jerry Lin 2023/12/06
%
%  data1:   datatable (require 'frame', no text column)
%  labelp:  list of markers (cell array)
%  k:       k nearst neighbor (by frame)    
%  sw1:     log switch

%% Initialization
tic;      
if nargin < 4
    sw1 = false;
end
data2 = varfun(@mean,data1,'GroupingVariables','frame');
data2.Properties.VariableNames = strrep(data2.Properties.VariableNames,'mean_','');
toc;

%% Caculate neighbors & testX & testY
if sw1
    [flag1,~] = knnsearch(log(data2{:,labelp}),log(data2{:,labelp}),'k',k);
else
    [flag1,~] = knnsearch(data2{:,labelp},data2{:,labelp},'k',k);
end
testX = mean(data2.Xt(flag1),2);
testY = mean(data2.Yt(flag1),2);
toc;

%% Output h1 (distance plot)

d1 = sqrt((testX-data2.Xt).^2+(testY-data2.Yt).^2);
figure,h2=scatter(data2.Xt,data2.Yt,10,log(d1),'fill');
%figure,h2=scatter(data2.Xt,data2.Yt,10,log(mean(d1,2)),'fill');
colormap(jet);colorbar;
set(gca,'XAxisLocation','bottom','YAxisLocation','left','ydir','reverse');
set(gca,'xticklabels',{});
set(gca,'yticklabels',{});
xlabel('X position');
ylabel('Y position');
daspect([1 1 1]);
set(gcf,'color','w');
title('Distance plot');
toc;

%% Output h2 (vector plot)
%dots = horzcat([data2.Xt, testX-data2.Xt],[data2.Yt,testY-data2.Yt]);

figure,scatter(data2.Xt,data2.Yt,5,1:length(data2.Xt));
hold on;
scatter(testX,testY,5,1:length(data2.Xt),'filled');
colormap(jet);
colorbar;
set(gca,'XAxisLocation','bottom','YAxisLocation','left','ydir','reverse');
set(gca,'xticklabels',{});
set(gca,'yticklabels',{});
daspect([1 1 1]);
set(gcf,'color','w');
xlabel('X position');
ylabel('Y position');
title('Vector plot');
toc;




