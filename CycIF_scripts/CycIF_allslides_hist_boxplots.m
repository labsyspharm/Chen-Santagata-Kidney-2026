

slideName = {'16868','23360','37932','54980'};

figure;
eval(strcat('temp1=sample',slideName{1},';'));
temp2 = NaN(length(temp1.X),length(slideName));

for slide = 1:length(slideName)
    myName = slideName{slide};
    eval(strcat('temp1=sample',myName,';'));
    data1 = temp1.CD3;
    data1 = asinh(data1(data1>0));
    %data1 = asinh(temp1.CD2-temp1.A488);
    subplot(2,2,slide);
    [pluscells, gate, peak,lowb,highb] = findgate2(data1,1,0.01);
    %[idx,c,gate, pluscells] = findgate2(data1,2,1);
    title(myName);
    temp2(1:length(data1),slide) = data1;
end
    
    
figure,boxplot(temp2);
set(gca,'xtick',1:length(slideName));
set(gca,'xticklabels',slideName);
set(gca,'XTickLabelRotation',90);
