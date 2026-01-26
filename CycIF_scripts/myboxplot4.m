function h1=myboxplot4(data1,group1,fontsize,grouporder)
%% My function boxplot3
%  Jerry Lin 2023/01/14
%  Usage myboxplot3(groupdata,groups,grouporder)

if nargin <=2
    fontsize=10;
end

if nargin <=3
    h1=boxplot(data1,group1);
else
    h1=boxplot(data1,group1,'GroupOrder',grouporder);
end    

if nargin <=3
    group4 = grp2idx(group1);
else
    group4 = zeros(length(group1));
    group4(ismember(group1,grouporder{1}))=1;
    group4(ismember(group1,grouporder{2}))=2;
    group4(ismember(group1,grouporder{3}))=3;
    group4(ismember(group1,grouporder{4}))=4;
end


ylims = ylim;
step1 = 0.1*(ylims(2)-ylims(1));
ylim([ylims(1),ylims(2)+4*step1]);
line([1,2],[ylims(2),ylims(2)],'LineWidth',0.5,'Color','k','Marker','|');
line([1,3],[ylims(2)+1.2*step1,ylims(2)+1.2*step1],'LineWidth',0.5,'Color','k','Marker','|');
line([1,4],[ylims(2)+2.2*step1,ylims(2)+2.2*step1],'LineWidth',0.5,'Color','k','Marker','|');

d1 = data1(group4==1);
d2 = data1(group4==2);
d3 = data1(group4==3);
d4 = data1(group4==4);

%-----p12------
[~,p12]=ttest2(d1,d2,'Vartype','unequal');
if p12 < 0.001
    clr1 = 'r';
    text(1.3,ylims(2)+0.5*step1,'p<0.001','FontSize',fontsize,'color',clr1);
else
    clr1 = 'k';
    text(1.3,ylims(2)+0.5*step1,strcat('p=',num2str(p12,'%0.3f')),'FontSize',fontsize,'color',clr1);
end
%-----p13------
[~,p13]=ttest2(d1,d3,'Vartype','unequal');
if p13 < 0.001
    clr1 = 'r';
    text(1.8,ylims(2)+1.8*step1,'p<0.001','FontSize',fontsize,'color',clr1);
else
    clr1 = 'k';
    text(1.8,ylims(2)+1.8*step1,strcat('p=',num2str(p13,'%0.3f')),'FontSize',fontsize,'color',clr1);
end
%-----p14------
[~,p14]=ttest2(d1,d4,'Vartype','unequal');
if p14 < 0.05
    clr1 = 'r';
    text(2.4,ylims(2)+2.8*step1,'p<0.001','FontSize',fontsize,'color',clr1);
else
    clr1 = 'k';
    text(2.4,ylims(2)+2.8*step1,strcat('p=',num2str(p14,'%0.3f')),'FontSize',fontsize,'color',clr1);
end

return;
