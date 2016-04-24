cd 'C:\Users\seth.koenig\Documents\MATLAB\SSCM\';
load('SSCM_ROI_DATA','ROI_All_xy','which_ROIs','replaced_ROIs','ROI_Cat_Labl','ROI_Cat');

attentionmaps = cell(90,4); %cols number of images, rows number of features
%feature 1: salience
%feature 2: monkeys
%feature 3: objects
%feature 4: red contrast

load('C:\Users\seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM\JN121212_2');
images1 = images;
fixationstats1 = fixationstats;
load('C:\Users\seth.koenig\Documents\MATLAB\IOR\Eye Data\SSCM\MP121206_2');
images2 = images;
fixationstats2 = fixationstats;
cd 'C:\Users\seth.koenig\Documents\MATLAB\SSCM\S2';

imageX = 800;
imageY = 600;
f = fspecial('gaussian',[256,256],24); %1 dva 2D gaussian filter
fixationpdf = cell(1,90);

for images = 1:90;
    reindexed = 2*images-1;
    fixationpdf{images} = zeros(600,800);
    
    fixations = fixationstats1{reindexed}.fixations;
    fixationtimes = fixationstats1{reindexed}.fixationtimes;
    if ~isempty(fixations)
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixations(:,1) = [];
        end
    end
    
    for fix = 1:size(fixations,2);
        fixx = round(fixations(1,fix));
        fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
        fixy = imageY - round(fixations(2,fix));
        fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
        fixationpdf{images}(fixy,fixx) = fixationpdf{images}(fixy,fixx)+1;
    end
    
    fixations = fixationstats2{reindexed}.fixations;
    fixationtimes = fixationstats2{reindexed}.fixationtimes;
    if ~isempty(fixations)
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixations(:,1) = [];
        end
    end
    
    for fix = 1:size(fixations,2);
        fixx = round(fixations(1,fix));
        fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
        fixy = imageY - round(fixations(2,fix));
        fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
        fixationpdf{images}(fixy,fixx) = fixationpdf{images}(fixy,fixx)+1;
    end
    
    fixationpdf{images} = imfilter(fixationpdf{images},f);
    fixationpdf{images} = fixationpdf{images}/sum(sum(fixationpdf{images}));
end
clear fixx fixy fix fixations fixationtimes

for image = 1:90;
    reindexed = 2*images-1;
    load([images1{reindexed} '-saliencemap'],'fullmap');
    fullmap = fullmap/sum(sum(fullmap));
    attentionmaps{image,1} = fullmap;
    img = imread([images1{reindexed} '.bmp']);
    Inorm = double(rgb2gray(img));
    Inorm(Inorm<0.1) = 1;
    r = double(img(:,:,1)); g = double(img(:,:,2)); b = double(img(:,:,3));
%     R = (r - (g+b)/2)./Inorm;
%     R = abs(R);
%     R = R-min(min((R))); 
% %     R = imfilter(R,f);
%     R = R/sum(sum(R));
    R = r./Inorm;
    R = R-min(min(R));
    R = imfilter(R,f);
    R = R/sum(sum(R));
    attentionmaps{image,4} = R;
end
clear r g b R Inorm fullmap

for images = 1:90;
    Objects = zeros(600,800);
    Monkeys = zeros(600,800);
    ROIS = which_ROIs{2}(images,:);
    replaced = replaced_ROIs{2}(images,:);
    for r = 1:length(ROIS);
        if ~isnan(ROIS(r));
            if r ~= replaced(2);
                   xy = ROI_All_xy{2}{ROIS(r)};
                   xy(:,2) = 600-xy(:,2);
                   ind = sub2ind([600,800],xy(:,2),xy(:,1));
                if ROI_Cat{2}{2}(ROIS(r),1) %if it is a monkey
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
    attentionmaps{images,2} = Monkeys;
    attentionmaps{images,3} = Objects;
end
clear Monkeys Objects ROIS replaced 
%% 
Bs = [];
for images = 1:90;
    y = fixationpdf{images}(1:end)';
    X = [];
    for feature = 1:size(attentionmaps,2);
       X = [X attentionmaps{image,feature}(1:end)'];
    end
    b = regress(y,X);
    Bs = [Bs b];
end
B_hat = mean(Bs,2);
B_hat = B_hat/sum(B_hat);
%%
X = [];
Y = [];
for images = 1:90;
    Y = [Y; fixationpdf{images}(1:end)'];
    x = [];
    for feature = 1:3%size(attentionmaps,2);
       x = [x attentionmaps{image,feature}(1:end)'];
    end
    X = [X;x];
end

%%
Y = [];
for images = 1:90;
   Y = [Y fixationpdf{images}(1:end)']; 
end
Y = Y(1:end);
%%
X = cell(1,4);
for images = 1:90;
    for c = 1:4
        X{c} = [X{c};attentionmaps{images,c}(1:end)'];
    end
end
X = cell2mat(X);
%%
B = glmfit(X,Y','normal','constant','off');
%% beta = mvregress(X,Y);