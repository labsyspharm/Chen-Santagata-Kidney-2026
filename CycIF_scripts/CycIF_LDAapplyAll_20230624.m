%% CycIF LDA apply all slides
%  Jerry Lin 2023/06/24
%
%  Require Variables: ldaAll
%  Need to run LDAcounts & LDAfitmodel first
%

%% Predict topics for all slides

tic;

for s=1:length(slideName)
    disp(strcat('Processing:',slideName{s}));
    data1 = eval(strcat('data',slideName{s}));
    
    
    count1 = varfun(@sum,data1,'GroupingVariables','frame','Inputvariables',labelp2);
    count1.topics = predict(ldaAll,count1{:,3:end});
    temp2 = count1(:,{'frame','topics'});
    if ismember('topics',data1.Properties.VariableNames)
        data1 = removevars(data1, 'topics');
    end
    data1 = join(data1,temp2,'keys','frame');
    
    toc;
    eval(strcat('data',slideName{s},'=data1;'));
end
