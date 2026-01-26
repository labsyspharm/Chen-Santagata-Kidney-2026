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
    data2 = data1(data1.ROIid>0,:);
    
    count1 = varfun(@sum,data2,'GroupingVariables','frame','Inputvariables',labelp2);
    count1.topics_ROI = predict(ldaROI,count1{:,3:end});
    temp2 = count1(:,{'frame','topics_ROI'});
    if ismember('topics_ROI',data1.Properties.VariableNames)
        data2 = removevars(data2, 'topics_ROI');
    end
    data2 = join(data2,temp2,'keys','frame');
    data1.topics_ROI = zeros(size(data1,1),1);
    data1.topics_ROI(ismember(data1.CellID,data2.CellID))=data2.topics_ROI;

    toc;
    eval(strcat('data',slideName{s},'=data1;'));
end
