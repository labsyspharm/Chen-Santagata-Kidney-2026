

figure;

for ch = 1:22;
        subplot(5,6,ch);
        data1 = alldata(:,ch);
        data1 = data1(1:5:end);
        %[idx,c]=kmeans(log(data1),2);
        histfit(log(data1+10),100,'kernel');
        title(channelnames(ch));
        %text(1,1,num2str(c));
        %xlim([4 10]);
    
end

