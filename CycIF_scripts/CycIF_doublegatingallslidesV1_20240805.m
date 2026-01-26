%% CycIF Doublegating all
%  Need slidename & doubeGates 
%  DoubleGating each slide independently
%  Jerry Lin  2024/08/05

%% -- Initialization --

doubleflag = true;
allmarkers = gateTable.Properties.VariableNames;



%% -- Gating all  & resampling --

for i =1:length(slideName)
        name1 = strcat('data',slideName{i});
        
        disp(strcat('Gating:',name1));
        data1 = eval(name1);
        disp(strcat('Processing:',slideName{i}));
        

        %-- Double gating --
        if doubleflag ==1
            for j=1:size(doubleGates,1)
                %gatename = strcat(doubleGates{i,1},'p',doubleGates{i,2},'p');
                gate1 = strcat(doubleGates{j,1},'p');
                gate2 = strcat(doubleGates{j,2},'p');
                gatename = strcat(gate1,gate2);
                data1{:,gatename}=data1{:,gate1} & data1{:,gate2};
            end
        end
        
        %-- re-assign --
        eval(strcat(name1,'=data1;'));
end

clear data1;


