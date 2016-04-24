function SSCM_SalienceIOR(FIXATIONFILE,img_dir,imageX,imageY)
% created by Seth Koenig 11/21/2012
% modified for SSCM content analysis and IOR 5/26/2014

% function determines rate of return fixations, the time between return
% fixations, time within trial of return, and salience at returned
% location. Inhibition of return was considered for mean fixation postions
% occruing within 0.5 dva of each other non-consequitively

% Inputs:
%   FIXATIONFILE: Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   ImageX,ImageY: x and y dimensions of images
%   Pairings: take all pairing or only closest unique pairing between first
%   and second fixation
%   img_dir: image director so can pull saliencemaps
%   imagefiles: image numbers so can pull salience maps
%   distancethreshold: distance categories for pairing fixations
%   novelconditions: novel condition numbers in fixationfile that represent novel
%   min_dist: minimum saccade distance out of area in dva   

%Outputs:
%   A .mat file named [FIXATIONFILE(1:end-13) '-SalienceIOR']
%   containg saccade statistics. See variable statvariablenames for
%   detailed explanation of variables in .mat file.

if nargin < 2
    error(['Not enough inputs: function requires FixationFile,'...
        'distance threhsold,imageX, imageY, and pairings.'])
end
if nargin < 3
    distancethreshold = [0;48];%in pixels 24 pixels/dva
end
if nargin < 5
    imageX = 800;
    imageY = 600;
end

load(FIXATIONFILE);


%fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal fixdist fix1dur
%fix2dur fixnum1 fixnum2 imagenum
min_dist = 10;
count = 1;
returnfixsal = NaN(500,15);

for cndlop=1:2:length(fixationstats)
    load([img_dir images{cndlop} '-saliencemap.mat'],'fullmap');
    imgnum = str2double(images{cndlop}(2:3));
    
    saliencemap = fullmap;
    fixations = fixationstats{cndlop}.fixations;
    if ~isempty(fixations)
        fixationtimes = fixationstats{cndlop}.fixationtimes;
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixationtimes(:,1) = [];
        end
        
        N=size(fixations,2);
        if N > 1
            [x,y]=meshgrid(1:N);
            i=find(ones(N)-eye(N)); %forms pairs except for self-pairing
            i=[x(i), y(i)];
            i(i(:,1) > i(:,2),:) = []; %repeat pairs
            i(i(:,1)+1 == i(:,2),:) = []; %removes consecutive in time pairs
            dist =sqrt((fixations(1,i(:,1))-fixations(1,i(:,2))).^2 +...
                (fixations(2,i(:,1))-fixations(2,i(:,2))).^2);
            wcount = 1;
            pairs = NaN(ceil(size(fixations,2)/2),3);
            while ~isempty(i);
                [minn,mind] = min(dist);
                minn = minn(1);
                mind = mind(1);
                
                middlefixes = i(mind,1)+1:i(mind,2)-1;
                middist = sqrt((fixations(1,i(mind,1))-fixations(1,middlefixes)).^2+...
                    (fixations(2,i(mind,2))-fixations(2,middlefixes)).^2);
                
                if any(middist >= min_dist*24)
                    tind = find(minn > distancethreshold(1,:) & minn <= distancethreshold(2,:));
                    if ~isempty(tind);
                        pairs(wcount,:) = [i(mind,:) tind];
                    end
                    [rmvind1,~] = find(i(:,1) == i(mind,1));
                    [rmvind2,~] = find(i(:,2) == i(mind,1));
                    [rmvind3,~] = find(i(:,1) == i(mind,2));
                    [rmvind4,~] = find(i(:,2) == i(mind,2));
                    rmvind = [rmvind1; rmvind2; rmvind3; rmvind4];
                    rmvind = unique(rmvind);
                    i(rmvind,:) = [];
                    dist(rmvind) = [];
                    wcount = wcount+1;
                else
                    dist(mind) = [];
                    i(mind,:) = [];
                end
            end
            pairs(isnan(pairs(:,1)),:) = [];
            
            for i = 1:size(pairs,1);
                
                spot = [ceil(fixations(:,pairs(i,1))) ceil(fixations(:,pairs(i,2)))];
                spot(2,:) = imageY-spot(2,:);
                spott = [(fixationtimes(1,pairs(i,1))+fixationtimes(2,pairs(i,1)))/2 ...
                    (fixationtimes(1,pairs(i,2))+fixationtimes(2,pairs(i,2)))/2];
                dist = sqrt((spot(1,1)-spot(1,2))^2+(spot(2,1)-spot(2,2))^2);
                spot(spot < 1) = 1;
                spot(1,spot(1,:) > imageX) = imageX;
                spot(2,spot(2,:) > imageY) = imageY;
                
                returnfixsal(count,:) = [...
                    spot(1,1) spot(2,1) spott(1) saliencemap(spot(2,1),spot(1,1))...
                    spot(1,2) spot(2,2) spott(2) saliencemap(spot(2,2),spot(1,2))...
                    dist diff(fixationtimes(:,pairs(i,1)))+1 ...
                    diff(fixationtimes(:,pairs(i,2)))+1 pairs(i,1) pairs(i,2) imgnum];
                %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal
                %fixdist fix1dur fix2dur fixnum1 fixnum2 imgnum

                count= count+1;
            end
        end
    end
end

returnfixsal(isnan(returnfixsal(:,1)),:) = [];

IORvariablenames = {
    'returnfixsal: [  %fix1x fix1y fix1t fix1sal fix2x fix2y fix2t fix2sal...';
    'fixdist fix1dur fix2dur fixnum1 fixnum2 imagenum]';
    };

save([FIXATIONFILE(1:10) '-SSCM_SalienceIOR.mat'],'returnfixsal','IORvariablenames','set')
end