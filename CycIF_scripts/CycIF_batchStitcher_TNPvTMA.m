%% Running Stitcher through a stitched plates to generate pyramids 
%  Jerry LIn 2019/12/27


%% Initialization

%allROI = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};

inputDIR = uigetdir('G:\TNP_2020\vTMA','Please choice the image directory');
outputDIR = uigetdir('G:\TNP_2020\vTMA\ometiff','Please choice the output directory');

ROIn = input('Please input # of images (10):');
tic;

for r=1:ROIn
    
    infile = strcat(inputDIR,'\TNP_vTMA',num2str(r),'.tif');
    disp(infile);
    outfile = strcat(outputDIR,'\TMA_vTMA',num2str(r),'.ome.tif');
    disp(outfile);
    myCommand = strcat({'"c:\Program Files\stitcher-0.4.2\bin\stitcher.bat" create_pyramid --input '},...
        infile,{' --output '},outfile);

    disp(myCommand);
    dos(myCommand{1},'-echo');
    toc;

end

