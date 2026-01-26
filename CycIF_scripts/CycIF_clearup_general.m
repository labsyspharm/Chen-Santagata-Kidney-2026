%% CycIF remove zeros, normalization & remove background
% Jerry LIn 2018/04/01
%
% require slideName variable from CycIF_import
% require function: CycIF_removezero, CycIF_normdata, CycIF_removebackgroud
%
%% Initialization

samplesize = input('Please input sample size:'); %5000;
flag1 = input('Normalization (Y/N)?','s');
flag2 = input('Convert to uint16 (Y/N)?','s');
flag3 = input('Remove background (Y/N)?','s');
sw1 = false;   % swtich for writing sample csv files

if samplesize >0
    allsample = table;
end

%% Remove zero from all data points

for i =1:length(slideName)
        name1 = strcat('data',slideName{i});
        disp(strcat('Removing zeros:',name1));
        data1 = eval(name1);
        data1 = CycIF_removezero(data1);
        eval(strcat(name1,'=data1;'));
end

clear data1;

%% X-axis normalization (flatten)

if ismember(flag1,{'Y','y'})
    for i =1:length(slideName)
            name1 = strcat('data',slideName{i});
            disp(strcat('Flatfield correction:',name1));
            data1 = eval(name1);   
            data1 = CycIF_normdata(data1);
            eval(strcat(name1,'=data1;'));
    end
    clear data1;
end

%% Convert to uint16

if ismember(flag2,{'Y','y'})
    for i =1:length(slideName)
            name1 = strcat('data',slideName{i});
            disp(strcat('Convert to uint16:',name1));
            data1 = eval(name1);   
            data1{:,2:end} = uint16(data1{:,2:end});
            eval(strcat(name1,'=data1;'));
    end
    clear data1;
end



%% Remove background & resampling

if ismember(flag3,{'y','Y'})
    for i =1:length(slideName)
            name1 = strcat('data',slideName{i});
            disp(strcat('Removing background & resampling:',name1));
            data1 = eval(name1);
            if ~exist('labels')
                labels = data1.Properties.VariableNames;
            end
            if ismember('A488',labels)
                data1 = CycIF_removebackground(data1,{'A488','A555'});
            else
                data1 = CycIF_removebackground(data1,{'FITC_1','Cy3_1'});
            end
            eval(strcat(name1,'=data1;'));
            % --- subsample ---
            if samplesize > 0
                sample1 = datasample(data1,samplesize);
                if(sw1)
                    filename = strcat('sample',slideName{i},'.csv');
                    writetable(sample1,filename);
                end

                eval(strcat('sample',slideName{i},'=sample1;'));
                sample1.slidename = repmat(slideName(i),length(sample1.X),1);
                if(isempty(allsample))
                    allsample = sample1;
                else
                    allsample = vertcat(allsample,sample1);
                end      
            end
            % ------------------
    end
    clear data1 sample1;
end
