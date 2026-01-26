%% CyCIF LDA model generation
%  jerry Lin 2023/06/24
%
%  Require variables: allcounts
%  Output:  lda1


%% --Calculate LDA & plots (markers only, All slides)-----

maxT = input('Please input number of topics (default=16):');
tic;
lda1 = fitlda(allcounts,maxT);
toc;
figure,imagesc(lda1.TopicWordProbabilities);colormap(jet);
title ('Topic Word probabilities');
set(gca,'ytick',1:length(labelp3));
set(gca,'yticklabels',labelp3);
xlabel('Topics');
colorbar;
caxis([0 0.25]);

figure
for topicIdx = 1:maxT
    subplot(4,ceil(maxT/4),topicIdx)
    temp1 = table;
    temp1.Word = labelp3;
    
    temp1.Count = lda1.TopicWordProbabilities(:,topicIdx);
    wordcloud(temp1,'Word','Count');
    title("Topic: " + topicIdx)
end
toc;
