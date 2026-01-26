function [h1,p]=myboxplot2(data1,group1,grouporder,fontsize)
%% My function boxplot2
%  Jerry Lin 2020/04/17
%  Usage myboxplot2(groupdata,groups,grouporder)
%  2021/09/19 add features to array only plot

if nargin <=1
    h1 = boxplot(data1);
elseif nargin <=2
    h1=boxplot(data1,group1);
else
    h1=boxplot(data1,group1,'GroupOrder',grouporder);
end
hold on;

ylims = ylim;
step1 = 0.1*(ylims(2)-ylims(1));
ylim([ylims(1),ylims(2)+2*step1]);
line([1,2],[ylims(2),ylims(2)],'LineWidth',1);

if nargin <=3
    fontsize = 9;
end

if nargin > 1
    group2 = grp2idx(group1);
    d1 = data1(group2==1);
    d2 = data1(group2==2);
else
    d1 = data1(:,1);
    d2 = data1(:,2);
end
% scatter(repmat(1,size(d1),1),d1,50,'k');
% scatter(repmat(2,size(d2),1),d2,50,'k');
hold off;
%[~,p]=ttest2(d1,d2,'Vartype','unequal');
[~,p]=ttest2(d1,d2,'Vartype','unequal');

if p<0.05
    color1 ='r';
    fontsize = fontsize+2;
else
    color1 = 'k';
end

if p<0.001
    text(1.25,ylims(2)+step1,'p<0.001','FontSize',fontsize,'Color',color1);
else
    text(1.25,ylims(2)+step1,strcat('p=',num2str(p,'%0.3f')),'FontSize',fontsize,'Color',color1);
end

return;
