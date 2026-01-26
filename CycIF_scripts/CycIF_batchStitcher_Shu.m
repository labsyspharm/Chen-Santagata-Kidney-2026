%% Running Stitcher through a stitched plates to generate pyramids 
%  Jerry LIn 2019/12/27


%% Initialization

rows = {'A','B','C','D','E','F','G','H'};
cols = {'01','02','03','04','05','06','07','08','09','10','11','12'};

inputDIR = uigetdir('H:\SHU\Raw_files\Cock','Please choice the image directory');
outputDIR = uigetdir('D:\Shu\OMETIFF\Cock','Please choice the output directory');

start_r = input('Please input starting row (1-8):');
end_r = input('Please input ending row (1-8):');
start_c = input('Please input starting column (1-12):');
end_c = input('Please input ending column (1-12):');
tic;

for r=start_r:end_r
    for c=start_c:end_c
    
    infile = strcat(inputDIR,'\',rows{r},cols{c},'.tif');
    disp(infile);
    outfile = strcat(outputDIR,'\',rows{r},cols{c},'.ome.tif');
    disp(outfile);
    myCommand = strcat({'"c:\Program Files\stitcher-0.4.2\bin\stitcher.bat" create_pyramid --input '},...
        infile,{' --output '},outfile);

    disp(myCommand);
    dos(myCommand{1},'-echo');
    toc;

    end
end

