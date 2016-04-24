% Code Imports tasks relevant data especially ROI data. Data is loaded and
% then processed into a more usable form, and saved for later use.
%
% WARNING this code currently guestimates the structure of the data. As I'm not
% the one that made this data and the code to analyze it was not available.
% I did my best to reasonably double check things, but it may not be perfect.
%
% Written by Seth Konig May 23, 2014

% [1] Import Data From Network
% [2] Organize ROIs by Image and Determine which ROIs are replaced
% [3] Minor Validation of the data above
% [4] Auxillary: Same as 2 but with plotting
% [5] Auxillary: Show Binary Matrix containing locations of ROI pixels
%% -----[1] Import Data From Network ---- %%
cd 'C:\Users\seth.koenig\Documents\MATLAB\SSCM';
dirdir= 'R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';

load([dirdir 'SSCM_Set1.2.3.4.5.6.7_indTrl.mat']);
load([dirdir 'SSCM_Set1.2.3.4.5.6.7_indTrl_ROI_Cat.mat']);
load([dirdir 'SSCM_Set1.2.3.4.5.6.7_ROI_All_xy.mat'])
load([dirdir 'SSCM_Set1.2.3.4.5.6.7_ROI_Cat.mat'])
load([dirdir 'SSCM_Set2.3.4.5.6.7_ROI_Face_xy.mat']);
load([dirdir 'SSCM_ROI_Scene.mat']);
load([dirdir 'SSCM_ROI_Include.mat']);
load('R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\SSCM_fruitInd_S1-7.mat');

clear indTrl indTrl indTrl_Labl include dirdir dateSaved time %I don't know what these are used for
cd 'C:\Users\seth.koenig\Documents\MATLAB\SSCM';
save('C:\Users\seth.koenig\Documents\MATLAB\SSCM\SSCM_ROI_DATA')
%% -----[2] Organize ROIs by Image and Determine which ROIs are replaced ---- %%
% assumes structure of ROIs that were replaced. Manually checked 100 or so
% images but it's possible that some ROIs were switched. At the very least
% the location of the ROIs (roughly their center) is conserved.
% 1 ROI was removed because of some bug that must have occured while it was
% created. This is the 13th ROI for the 84th image in Set3. Other irregularities
% may be present Therefore therefore I removed this image from analysis.
cd 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM';
imgdir = 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM\';

load('SSCM_ROI_DATA')
replaced_ROIs = cell(1,7); %which ROIs are associated with replaced objects
which_ROIs = cell(1,7);  %which ROIs go with which image
for Set = 2:7;
    cd([imgdir 'S' num2str(Set) '\']);
    files = dir('*.bmp'); %find all the file names
    files = {files.name}';
    
    replaced_ROIs{Set} = NaN(size(which_ROIs,1),2); %which ROIs are associated with replaced objects
    imagenumber = NaN(1,length(files));
    for f = 1:length(files);
        imagenumber(f) = str2double(files{f}(4:5));
    end
    
    which_ROIs{Set} = NaN(90,13); %which ROIs go with which image
    for im = 1:90;
        ind = find(ROI_Scene{Set} == im);
        which_ROIs{Set}(im,1:length(ind)) = ind;
    end
    if Set == 3
        which_ROIs{Set}(84,:) = NaN; %at least 1 weird ROI so remove all
    end
    
    %try to find overlapping ROIs to determine which ones were
    %replaced since their location should be approximately the same
    for im = 1:max(imagenumber);
        imagefiles_num = find(imagenumber == im);
        if length(imagefiles_num ) == 2; %if only 1 image, else have replaced trials if there are 2
            if sum(~isnan(which_ROIs{Set}(im,:))) == 13 %double check since must be a replaced image Set
                Rnum = 1;
                not_found = true;
                while Rnum <= 12 && not_found
                    %make binary matrix with 1's where the ROI is and 0's where it's not
                    temp1 = NaN(600,800);
                    xy = ROI_All_xy{Set}{which_ROIs{Set}(im,Rnum)};
                    xy(xy == 0) = 1;
                    temp1(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
                    for comparisonROI = Rnum+1:13
                        %make another binary matrix with 1's where ROI is for the ROI we're comparing to see if ROI1 overalps with ROI2
                        temp2 = NaN(600,800);
                        xy = ROI_All_xy{Set}{which_ROIs{Set}(im,comparisonROI)};
                        xy(xy == 0) = 1;
                        temp2(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
                        if any(any(temp1 == temp2)) %if they overlap at all, NaN ~= NaN!
                            minsize = min(nansum(nansum(temp1)),nansum(nansum(temp2)));
                            if sum(sum(temp1 == temp2)) > 0.33*minsize %to ensure ROIs that just touch don't get falsely flagged
                                %ordering appears to always be the same as far as I can tell
                                replaced_ROIs{Set}(im,2) = Rnum; %object to be replaced
                                replaced_ROIs{Set}(im,1) = comparisonROI; %object that which replaced
                                not_found = false; %to break while loop
                                break; %to break for loop
                            end
                        end
                    end
                    Rnum=Rnum+1;
                end
            end
        end
    end
end
save('C:\Users\seth.koenig\Documents\MATLAB\SSCM\SSCM_ROI_DATA','replaced_ROIs','which_ROIs','-append')
%% -----[3] Minor Validation of the data above---- %%
% not the best validation but a good check nonetheless. Determines if there
% are 60 replaced objects detected per set as there should be. Also, it
% checks if the 8th and 9th ROI are correctly associated with the replaced
% object as it seems these the only ROIs that were changed.
number_replaced = [];
for Set = 2:7
    number_replaced(Set) = sum(~isnan(replaced_ROIs{Set}(:,1)));
end
number_replaced(1) = [];
number_replaced(2) = number_replaced(2)+1; %removed 84th image from set 3 on purpose
if any(number_replaced ~= 60)
    disp('Error Likely Occured in detecting replaced ROI. Double check answer')
end

for Set = 2:7
    rr = replaced_ROIs{Set};
    rr(isnan(rr)) = [];
    %figure
    %hist(rr) %can display if all images are 8 or 9 and should be uniform
    if any(rr ~= 8 & rr ~= 9)
        oddind = find(rr ~= 8 & rr ~= 9);
        for odd = 1:length(oddind);
            [img,ROIind] = find(replaced_ROIs{Set} == rr(oddind(odd)));
            disp('Odd ROI detected as being replaced');
            disp(['Set #' num2str(Set)])
            disp(['Image #: ' num2str(img) ' ROI # ' num2str(rr(oddind(odd)))])
        end
    end
end
%% -----[4] Auxillary: Same as 2 but with plotting ---- %%
% can be used to visually double check if the ROIs are correclty
% interpreted and that you have correctly organized and processed them. 
cd 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM';
imgdir = 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM\';

clr = ['r*';'g*';'c*';'m*'];

load('SSCM_ROI_DATA')
replaced_ROIs = cell(1,7); %which ROIs are associated with replaced objects
which_ROIs = cell(1,7);  %which ROIs go with which image
for Set = 2%2:7;
    cd([imgdir 'S' num2str(Set) '\']);
    files = dir('*.bmp'); %find all the file names
    files = {files.name}';
    
    replaced_ROIs{Set} = NaN(size(which_ROIs,1),2); %which ROIs are associated with replaced objects
    imagenumber = NaN(1,length(files));
    for f = 1:length(files);
        imagenumber(f) = str2double(files{f}(4:5));
    end
    
    which_ROIs{Set} = NaN(90,13); %which ROIs go with which image
    for im = 1:90;
        ind = find(ROI_Scene{Set} == im);
        which_ROIs{Set}(im,1:length(ind)) = ind;
    end
    if Set == 3
        which_ROIs{Set}(84,:) = NaN; %at least 1 weird ROI so remove all
    end
    
    %try to find overlapping ROIs to determine which ones were
    %replaced since their location should be approximately the same
    
    %replaced_ROI_ordering = NaN(1,size(which_ROIs,1)); %determine which of the replaced ROIs is first
    for im = 1:5%:max(imagenumber);
        imagefiles_num = find(imagenumber == im);
        if length(imagefiles_num ) == 2; %if only 1 image, else have replaced trials if there are 2
            figure
            img1 = imread(files{imagefiles_num(1)});
            
            subplot(1,2,1)
            image(img1);
            axis off
            
            img2 = imread(files{imagefiles_num(2)});
            subplot(1,2,2)
            image(img2);
            axis off
            
            if sum(~isnan(which_ROIs{Set}(im,:))) == 13 %double check since must be a replaced image Set
                Rnum = 1;
                not_found = true;
                while Rnum <= 12 && not_found
                    %make binary matrix with 1's where the ROI is and 0's where it's not
                    temp1 = NaN(600,800);
                    xy = ROI_All_xy{Set}{which_ROIs{Set}(im,Rnum)};
                    xy(xy == 0) = 1;
                    temp1(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
                    for comparisonROI = Rnum+1:13
                        %make another binary matrix with 1's where ROI is for the ROI we're comparing to see if ROI1 overalps with ROI2
                        temp2 = NaN(600,800);
                        xy = ROI_All_xy{Set}{which_ROIs{Set}(im,comparisonROI)};
                        xy(xy == 0) = 1;
                        temp2(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
                        if any(any(temp1 == temp2)) %if they overlap at all, NaN ~= NaN!
                            minsize = min(nansum(nansum(temp1)),nansum(nansum(temp2)));
                            if sum(sum(temp1 == temp2)) > 0.33*minsize %to ensure ROIs that just touch don't get falsely flagged
                                %ordering appears to always be the same as far as I can tell
                                replaced_ROIs{Set}(im,2) = Rnum; %object to be replaced
                                replaced_ROIs{Set}(im,1) = comparisonROI; %object that which replaced
                                not_found = false; %to break while loop
                                break; %to break for loop
                            end
                        end
                    end
                    Rnum=Rnum+1;
                end
            end
            
             for imgnum = 1:2
                subplot(1,2,imgnum)
                hold on
                for r = 1:size(which_ROIs{Set},2);
                    ROInum = which_ROIs{Set}(im,r);
                    xy = ROI_All_xy{Set}{ROInum};
                    xy(xy == 0) = 1;
                    temp = zeros(600,800);
                    temp(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
                    temp = im2bw(temp);
                    temp = imfill(temp,'holes');
                    border_points = Find_Edges(temp);
                    if ROI_Cat{Set}{2}(r,1) == 1;
                        type = clr(1,:); %monkeys
                    else
                        type = clr(2,:);%objects
                    end
                    if fruitInd{Set}(which_ROIs{Set}(im,r)) == 1;
                        type = clr(3,:);%fruit objects
                    end
                    if any(replaced_ROIs{Set}(im,:) == r);
                        if find(replaced_ROIs{Set}(im,:) == r) == imgnum;
                            type =clr(4,:); %for replaced objects
                        else
                            continue %skip this ROI and go to the next iterion or just end
                        end
                    end
                    plot(border_points(:,2),600-border_points(:,1),type)
                end
             end
            
             % Can be used to manually switch which ROI was replaced if you
             % get it wrong. Note: Manually swtiching was not previously
             % implemented. 
             %str = input('Are the replaced ROIs in the Correct order?','s');
             %if strcmpi(str,'y')
             %  replaced_ROI_ordering(im) = true;
             %else
             %  replaced_ROI_ordering(im) = false;
             %end
            
        end
    end
end


%% -----[5] Auxillary: Show Binary Matrix containing locations of ROI pixels ---- %%
%similar to the code above but shows a matrix containing 1's where ROIs are
%located and 0's where there is background. 

%pic and set and image to explore
Set = 3;
im = 84; %image number

cd([imgdir 'S' num2str(Set) '\']);
files = dir('*.bmp');
files = {files.name}';

replaced_ROIs{Set} = NaN(size(which_ROIs,1),2); %which ROIs are associated with replaced objects

imagenumber = NaN(1,length(files));
for f = 1:length(files);
    imagenumber(f) = str2double(files{f}(4:5));
end

which_ROIs = NaN(90,13); %which ROIs go with which image
for im = 1:90;
    ind = find(ROI_Scene{Set} == im);
    which_ROIs(im,1:length(ind)) = ind;
end

imagefiles_num = find(imagenumber == im);
if length(imagefiles_num ) == 2; %if only 1 image, else have replaced trials if there are 2
    %try to find overlapping ROIs to determine which ones were
    %replaced since their location should be approximately the same
    if sum(~isnan(which_ROIs(im,:))) == 13 %double check since must be a replaced image Set
        Rnum = 1;
        not_found = true;
        while Rnum <= 12 && not_found
            %make binary matrix with 1's where the ROI is and 0's where it's not
            temp1 = NaN(600,800);
            xy = ROI_All_xy{Set}{which_ROIs(im,Rnum)};
            xy(xy == 0) = 1;
            temp1(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
            for comparisonROI = Rnum+1:13
                %make another binary matrix with 1's where ROI is for the ROI we're comparing to see if ROI1 overalps with ROI2
                temp2 = NaN(600,800);
                xy = ROI_All_xy{Set}{which_ROIs(im,comparisonROI)};
                xy(xy == 0) = 1;
                temp2(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
                if any(any(temp1 == temp2)) %if they overlap at all, NaN ~= NaN!
                    minsize = min(nansum(nansum(temp1)),nansum(nansum(temp2)));
                    if sum(sum(temp1 == temp2)) > 0.33*minsize %to ensure ROIs that just touch don't get falsely flagged
                        %ordering appears to always be the same as far as I can tell
                        replaced_ROIs{Set}(im,2) = Rnum; %object to be replaced
                        replaced_ROIs{Set}(im,1) = comparisonROI; %object that which replaced
                        not_found = false; %to break while loop
                        break; %to break for loop
                    end
                end
            end
            Rnum=Rnum+1;
        end
    end
    
    for imgnum = 1:2
        subplot(1,2,imgnum)
        temp = zeros(600,800);
        for r = 1:size(which_ROIs,2)-1;
            ROInum = which_ROIs(im,r);
            xy = ROI_All_xy{Set}{ROInum};
            xy(xy == 0) = 1;
            temp(sub2ind([600,800],xy(:,2),xy(:,1))) = 1;
        end
        imagesc(temp)
    end
end