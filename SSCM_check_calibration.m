% eye data was stored else where because I originally used it for another analysis
% currently data is for saline delivered trials not Oxytocin trials
eye_data_dir = 'C:\Users\Seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM\';
SSCM_eye_dir = 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM\EyeTraces\';
image_dir = 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM\';

eyedat_files = {'JN121212_2','MP121206_2';
    'JN121115_2','MP130705_2';
    'JN130514_2','MP121210_2';
    'JN121119_2','MP130709_2';
    'JN130516_2','MP121212_2';
    'JN121121_2','MP130711_2'};

imageX = 800; %horizontal size of images
imageY = 600; %vertical size of images
for set = 2:7;
    cd([image_dir 'S' num2str(set)])
    mkdir([SSCM_eye_dir 'S' num2str(set) '\JN\']);
    mkdir([SSCM_eye_dir 'S' num2str(set) '\MP\']);
    
    load([eye_data_dir eyedat_files{set-1,1}]);
    fixationstats1 = fixationstats;
    for imgnum = 1:2:length(fixationstats1)
        img = imread([images{imgnum} '.bmp']);
        figure
        imshow(img);
        hold on
        plot(fixationstats1{imgnum}.XY(1,:),imageY-fixationstats1{imgnum}.XY(2,:),'r')
        hold off
        axis off
        print('-djpeg',[SSCM_eye_dir 'S' num2str(set) '\JN\' images{imgnum}],'-r100');
        close
    end
    
    load([eye_data_dir eyedat_files{set-1,2}]);
    fixationstats2 = fixationstats;
    for imgnum = 1:2:length(fixationstats2);
        img = imread([images{imgnum} '.bmp']);
        figure
        imshow(img);
        hold on
        plot(fixationstats2{imgnum}.XY(1,:),imageY-fixationstats2{imgnum}.XY(2,:),'r')
        hold off
        axis off
        print('-djpeg',[SSCM_eye_dir 'S' num2str(set) '\MP\' images{imgnum}],'-r100');
        close
    end
end



