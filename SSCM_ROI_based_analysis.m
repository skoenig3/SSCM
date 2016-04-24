% eye data was stored else where because I originally used it for another analysis
% currently data is for saline delivered trials not Oxytocin trials
eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM\';
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SSCM\';

cd(image_dir)
load('SSCM_ROI_DATA','ROI_All_xy','which_ROIs','replaced_ROIs','ROI_Cat_Labl','ROI_Cat');

cross_feature_ROCs = cell(6,90,6); %set number x image number x feature number
%feature 1: salience
%feature 2: monkeys
%feature 3: objects
%feature 4: red contrast
%feature 5: all features combined
%feature 6: semantic features only (i.e. no salience)

B_hat = [0.6355; 0.1765; 0.1404; 0.0476]; %briefly guessed from means of linear regression on S2
% B_hat = [0.3 0.3 0.1 0.3]; 

eyedat_files = {'JN121212_2','MP121206_2';
    'JN121115_2','MP130705_2';
    'JN130514_2','MP121210_2';
    'JN121119_2','MP130709_2';
    'JN130516_2','MP121212_2';
    'JN121121_2','MP130711_2'};

imageX = 800; %horizontal size of images
imageY = 600; %vertical size of images
f = fspecial('gaussian',[256,256],24); %1 dva 2D gaussian filter

clrs = ['rgbmck'];
% figure
% hold on

for set = 2:7;
    cd([image_dir 'S' num2str(set)])
    
    load([eye_data_dir eyedat_files{set-1,1}]);
    fixationstats1 = fixationstats;
%     load([eye_data_dir eyedat_files{set-1,2}]);
%     fixationstats2 = fixationstats;
%     
    for imgnum = 1:floor(min(length(fixationstats1)/2));%,length(fixationstats2))/2);
        disp(['set #' num2str(set) ' Image #' num2str(imgnum)])
        reindexed = 2*imgnum-1; %reindexed cuz every other is repeat/replaced trial
        
        %fullmap contains salience map for that image
        load([images{reindexed} '-saliencemap'],'fullmap');
        fullmap = fullmap/sum(sum(fullmap));
        
        
        img = imread([images{reindexed} '.bmp']);
        Inorm = double(rgb2gray(img)); %gray scale image intensity
        Inorm(Inorm<0.1) = 1; %don't normalize by low image intensity
        r = double(img(:,:,1)); %red layer of image
        %         g = double(img(:,:,2)); %green layer of image
        %         b = double(img(:,:,3)); %blue layer of image
        %         R = (r - (g+b)/2)./Inorm; %normalized red hue
        %
        %         %turn red hue into pdf
        %         R = abs(R);
        %         R = R-min(min((R)));
        %         R = R/sum(sum(R));
        R = r./Inorm;
        R = R-min(min(R));
        R = imfilter(R,f);
        R = R/sum(sum(R));
        
        Objects = zeros(600,800);%will hold location of monkeys
        Monkeys = zeros(600,800);
        ROIS = which_ROIs{set}(imgnum,:); %the ordinal number of the ROI
        replaced = replaced_ROIs{set}(imgnum,:); %col1 is which ROI is to be replaced, col2 is the roi that replaced
        for r = 1:length(ROIS);
            if ~isnan(ROIS(r));
                if r ~= replaced(2);
                    xy = ROI_All_xy{set}{ROIS(r)}; %pixel cordinates of the roi
                    xy(:,2) = 600-xy(:,2); %flip vertically to align with everything
                    ind = sub2ind([600,800],xy(:,2),xy(:,1));
                    if ROI_Cat{set}{2}(ROIS(r),1) %if it is a monkey
                        Monkeys(ind) = 1;
                    else
                        Objects(ind) = 1;
                    end
                end
            end
        end
        Monkeys = imfilter(Monkeys,f);
        Objects = imfilter(Objects,f);   
        Monkeys = Monkeys/sum(sum(Monkeys));
        Objects = Objects/sum(sum(Objects));
                
        combined_features = B_hat(1)*fullmap+B_hat(2)*Monkeys+B_hat(3)*Objects+B_hat(4)*R;%linear weighted sum
        combined_semantic = Monkeys+Objects; %weighting currently unknown
        
        for feature =1:6;
            if feature == 1
                temp = fullmap/max(max(fullmap));
            elseif feature == 2
                temp = Monkeys/max(max(Monkeys));
            elseif feature == 3
                temp = Objects/max(max(Objects));
            elseif feature == 4
                temp = R/max(max(R));
            elseif feature == 5
                temp = combined_features/max(max(combined_features));
            elseif feature == 6
                temp = combined_semantic/max(max(combined_semantic));
            end
            
            fixations = fixationstats1{reindexed}.fixations; %fixation locations from monkey 1
            %remove first fixation if it is the crosshair fixation
            if ~isempty(fixations)
                if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                        fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
                    fixations(:,1) = [];
                    fixations(:,1) = [];
                end
            end
            
            %calculate feature values at observed fixation locations and
            %random "fixation" locations for monkey 1
            values = [];
            random_values = [];
            for fix = 1:size(fixations,2);
                fixx = round(fixations(1,fix));
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - round(fixations(2,fix));
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                
                ry = randi(imageY); ry(ry < 1) = 1;
                rx = randi(imageX); rx(rx < 1) = 1;
                
                values = [values temp(fixy,fixx)];
                random_values = [random_values temp(ry,rx)];
            end
            
%             fixations = fixationstats2{reindexed}.fixations;%fixation locations from monkey 2
%             %remove first fixation if it is the crosshair fixation
%             if ~isempty(fixations)
%                 if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
%                         fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
%                     fixations(:,1) = [];
%                     fixations(:,1) = [];
%                 end
%             end
%             
%             %calculate feature values at observed fixation locations and
%             %random "fixation" locations for monkey 2
%             for fix = 1:size(fixations,2);
%                 fixx = round(fixations(1,fix));
%                 fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
%                 fixy = imageY - round(fixations(2,fix));
%                 fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
%                 
%                 ry = randi(imageY); ry(ry < 1) = 1;
%                 rx = randi(imageX); rx(rx < 1) = 1;
%                 
%                 values = [values temp(fixy,fixx)];
%                 random_values = [random_values temp(ry,rx)];
%             end
            
            len = length(values);
            thresh = 0:0.01:1;
            TP = NaN(1,length(thresh)); %True positive
            FA = NaN(1,length(thresh)); %False alarm
            for ii = 1:length(thresh)
                TP(ii) = sum(values >= thresh(ii))/len;
                FA(ii) = sum(random_values >= thresh(ii))/len;
            end
            cross_feature_ROCs{set-1,imgnum,feature} = [TP;FA];
            %             plot(FA,TP,clrs(feature))
            %             pause(0.05)
        end
    end
end
cd('C:\Users\seth.koenig\Documents\MATLAB\SSCM\')
save('ROCs.mat','cross_feature_ROCs');
%%
cd('C:\Users\seth.koenig\Documents\MATLAB\SSCM')
load('ROCs.mat')


AUROC = NaN(6,90,6);
all_ROC = cell(2,6);
for set = 2:7
    for imgnum = 1:90;
        for features = 1:6
            if ~isempty(cross_feature_ROCs{set-1,imgnum,features})
                TP = cross_feature_ROCs{set-1,imgnum,features}(1,:);
                FA = cross_feature_ROCs{set-1,imgnum,features}(2,:);
                all_ROC{1,features} = [all_ROC{1,features}; TP];
                all_ROC{2,features} = [all_ROC{2,features}; FA];
                AUROC(set-1,imgnum,features) = -trapz(FA,TP);
            end
        end
    end
end

clrs = ['rgbmck'];
figure
hold on
for features = 1:6;
    avg_TP = mean(all_ROC{1,features});
    avg_FA = mean(all_ROC{2,features});
    plot(avg_FA,avg_TP,clrs(features))
end
plot([0 1],[0 1],'k--','linewidth',3)%unity line
axis square
xlabel('False Alarm Rate')
ylabel('True Positive Rate')
legend('Salience-only','Monkeys-only','Objects-only','Red-only','All Combined','Monkeys+Objects');

A = [];
A = reshape(AUROC,[6*90,6]);
nanmean(A,1) %average ROC values

%%
% eye data was stored else where because I originally used it for another analysis
% currently data is for saline delivered trials not Oxytocin trials
eye_data_dir = 'C:\Users\Seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM\';
image_dir = 'C:\Users\Seth.koenig\Documents\MATLAB\SSCM\';

cd(image_dir)
load('SSCM_ROI_DATA','ROI_All_xy','which_ROIs','replaced_ROIs','ROI_Cat_Labl','ROI_Cat','fruitInd');


eyedat_files = {'JN121212_2','MP121206_2';
    'JN121115_2','MP130705_2';
    'JN130514_2','MP121210_2';
    'JN121119_2','MP130709_2';
    'JN130516_2','MP121212_2';
    'JN121121_2','MP130711_2'};

imageX = 800; %horizontal size of images
imageY = 600; %vertical size of images

% Fixation Location Content Code
% 0: background
% 1: Monkey
% 2: Object
% 3: Fruit

Return_Location_Content = cell(6,2,3); %Set by monkey by ROI extension size (0,12 pixels, 24 pixels)
% ROI extension may account for any calibration issues
for Set = 2:7
    for monk = 1:2
        for ext = 1:3
            Return_Location_Content{Set-1,monk,ext} = NaN(500,3);%content by fixation number and by it's salience
        end
    end
end
count = ones(6,2,3);
for Set = 2:7;
    
    %for JN
    load([eye_data_dir eyedat_files{Set-1,1} '-SSCM_SalienceIOR']);
    returnfixsal1 = returnfixsal;
    %for MP
    load([eye_data_dir eyedat_files{Set-1,2} '-SSCM_SalienceIOR']);
    returnfixsal2 = returnfixsal;
    
    %create matrices labeling where objects and monkeys are  in a scene
    for imgnum = 1:90
        disp(['Set #' num2str(Set) ' Image #' num2str(imgnum)])
        Content = zeros(600,800);%will hold location of monkeys all content
        Content_Extended1 = zeros(600,800); %content extended by 24 pixels/1 dva
        Content_Extended2 = zeros(600,800); %content extended by 48 pixels/2 dva
        ROIS = which_ROIs{Set}(imgnum,:); %the ordinal number of the ROI
        replaced = replaced_ROIs{Set}(imgnum,:); %col1 is which ROI is to be replaced, col2 is the roi that replaced
        for r = 1:length(ROIS);
            if ~isnan(ROIS(r));
                if r ~= replaced(2);
                    xy = ROI_All_xy{Set}{ROIS(r)}; %pixel cordinates of the roi
                    xy(:,2) = imageY-xy(:,2); %flip vertically to align with everything
                    xy(xy < 1) = 1;
                    ind = sub2ind([600,800],xy(:,2),xy(:,1));
                    [extind1] = extend_ROI(xy,24,imageX,imageY);
                    [extind2] = extend_ROI(xy,48,imageX,imageY);
                    if ROI_Cat{Set}{2}(ROIS(r),1) %if it is a monkey
                        Content(ind) = 1;
                        Content_Extended1(extind1) = 1;
                        Content_Extended2(extind2) = 1;
                    else
                        if fruitInd{Set}(ROIS(r))
                            Content(ind) = 3;
                            Content_Extended1(extind1) = 3;
                            Content_Extended2(extind2) = 3;
                        else
                            Content(ind) = 2;
                            Content_Extended1(extind1) = 2;
                            Content_Extended2(extind2) = 2;
                        end
                    end
                end
            end
        end
        
        %For JN
        for_this_img = find(returnfixsal1(:,end) == imgnum);
        if ~isempty(for_this_img);
            for f = 1:length(for_this_img);
                %return locations
                fixx = returnfixsal1(for_this_img(f),5);
                fixy = returnfixsal1(for_this_img(f),6);
                fixnum = returnfixsal1(for_this_img(f),12);
                sal = returnfixsal1(for_this_img(f),8);
                Return_Location_Content{Set-1,1,1}(count(Set-1,1,1),1) = Content(fixy,fixx);
                Return_Location_Content{Set-1,1,2}(count(Set-1,1,2),1) = Content_Extended1(fixy,fixx);
                Return_Location_Content{Set-1,1,3}(count(Set-1,1,3),1) = Content_Extended2(fixy,fixx);
                Return_Location_Content{Set-1,1,1}(count(Set-1,1,1),2) = fixnum;
                Return_Location_Content{Set-1,1,2}(count(Set-1,1,2),2) = fixnum;
                Return_Location_Content{Set-1,1,3}(count(Set-1,1,3),2) = fixnum;
                Return_Location_Content{Set-1,1,1}(count(Set-1,1,1),3) = sal;
                Return_Location_Content{Set-1,1,2}(count(Set-1,1,2),3) = sal;
                Return_Location_Content{Set-1,1,3}(count(Set-1,1,3),3) = sal;
                count(Set-1,1,1) = count(Set-1,1,1)+1;
                count(Set-1,1,2) = count(Set-1,1,2)+1;
                count(Set-1,1,3) = count(Set-1,1,3)+1;
            end
        end
        
        %For MP
        for_this_img = find(returnfixsal2(:,end) == imgnum);
        if ~isempty(for_this_img);
            for f = 1:length(for_this_img);
                %return locations
                fixx = returnfixsal2(for_this_img(f),5);
                fixy = returnfixsal2(for_this_img(f),6);
                fixnum = returnfixsal2(for_this_img(f),12);
                sal = returnfixsal2(for_this_img(f),8);
                Return_Location_Content{Set-1,2,1}(count(Set-1,2,1),1) = Content(fixy,fixx);
                Return_Location_Content{Set-1,2,2}(count(Set-1,2,2),1) = Content_Extended1(fixy,fixx);
                Return_Location_Content{Set-1,2,3}(count(Set-1,2,3),1) = Content_Extended2(fixy,fixx);
                Return_Location_Content{Set-1,2,1}(count(Set-1,2,1),2) = fixnum;
                Return_Location_Content{Set-1,2,2}(count(Set-1,2,2),2) = fixnum;
                Return_Location_Content{Set-1,2,3}(count(Set-1,2,3),2) = fixnum;
                Return_Location_Content{Set-1,2,1}(count(Set-1,2,1),3) = sal;
                Return_Location_Content{Set-1,2,2}(count(Set-1,2,2),3) = sal;
                Return_Location_Content{Set-1,2,3}(count(Set-1,2,3),3) = sal;
                count(Set-1,2,1) = count(Set-1,2,1)+1;
                count(Set-1,2,2) = count(Set-1,2,2)+1;
                count(Set-1,2,3) = count(Set-1,2,3)+1;
            end
        end
        
    end
end
%%
save('ReturnData','Return_Location_Content')
content1=cell(1,3);
content2 = cell(1,3);

for Set = 2:7
    for monk = 1:2
        for ext = 1:3
            if monk == 1;
                content1{ext} = [content1{ext}; Return_Location_Content{Set-1,monk,ext}];
            else
                content2{ext} = [content1{ext}; Return_Location_Content{Set-1,monk,ext}];
            end
        end
    end
end

figure
subplot(1,3,1)
hist(content1{1})
title('No extension')
subplot(1,3,2)
hist(content1{2})
title('1/2 dva extension')
subplot(1,3,3)
hist(content1{3})
title('1 dva extension')
subtitle('JN')


figure
subplot(1,3,1)
hist(content2{1})
title('No extension')
subplot(1,3,2)
hist(content2{2})
title('1/2 dva extension')
subplot(1,3,3)
hist(content2{3})
title('1 dva extension')
subtitle('MP')
%%
%get percent of returns

figure
subplot(1,2,1)
hold on
bar(1,length(find(content1{1}(:,1) == 0))/length(content1{1}),'r')
bar(2,length(find(content1{1}(:,1) == 1))/length(content1{1}),'g')
bar(3,length(find(content1{1}(:,1) == 2))/length(content1{1}),'b')
bar(4,length(find(content1{1}(:,1) == 3))/length(content1{1}),'m')
title('JN')
ylabel('Percent Returns');
set(gca,'Xtick',1:4);
set(gca,'XtickLabel',{'"Background"','Monkeys','Objects','Fruit'})

subplot(1,2,2)
hold on
bar(1,length(find(content2{1}(:,1) == 0))/length(content2{1}),'r')
bar(2,length(find(content2{1}(:,1) == 1))/length(content2{1}),'g')
bar(3,length(find(content2{1}(:,1) == 2))/length(content2{1}),'b')
bar(4,length(find(content2{1}(:,1) == 3))/length(content2{1}),'m')
title('MP')
ylabel('Percent Returns')
set(gca,'Xtick',[1 2 3 4]);
set(gca,'XtickLabel',{'"Background"','Monkeys','Objects','Fruit'})
%%
%how salience at the returns to background or other content
backgroundind1 = find(content1{1}(:,1) == 3);
backgroundind2 = find(content2{2}(:,1) == 3);

sal1 = content1{1}(backgroundind1,3);
sal2 = content2{2}(backgroundind2,3);
nanmean(sal1)
nanmean(sal2)