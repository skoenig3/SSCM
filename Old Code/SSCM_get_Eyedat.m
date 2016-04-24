%% Load Behavioral Data
% vardir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/';
% inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';
dirRoot='C:\Users\drew\Documents\Buffalo Rotation\Scene Manipulation\';
dirDat=[dirRoot 'Data Files\'];
dirItm=[dirRoot 'Itm Files\'];
dirCnd=[dirRoot 'Cnd Files\'];
dirFig=[dirRoot 'Figures\'];
dirVar=[dirRoot 'Variables\'];
dirInd=[dirVar 'ROI_Ind\Files\'];
dirEye=[dirVar 'ROI_Ind\Files\eyedat\'];
dirFix=[dirVar 'ROI_Ind\Files\datFix\'];
dirSave=[dirFig 'Rank Images\'];  

dateSaved=date;
time=clock;
time=[num2str(time(4)),'h',num2str(time(5)),'m'];

expectedCnd=1010:1189;

% Number of Trials?
T=180;
% What subjects have seen each set?
sList={...
    %Set1
    {'JN','MP','IW','TT'};...
    %Set2
    {'JN','MP'};...
    %Set3
    {'JN','MP'}...
    %Set4
    {'JN','MP'}...
    %Set5
    {'JN','MP'}...
    %Set6
    {'JN','MP'}...
    %Set7
    {'JN','MP'}...
    };

% Data Files separated by set and subject
fList={...
    %Set1
    {{'JN120717.4','JN130513.3'},{'MP120716.1','MP130703.4'},{'IW120425.1'},{'TT120717.5'}}...
    %Set2
    {{'JN121108.2','JN121212.2'},{'MP121206.2','MP130704.2'}}...
    %Set3
    {{'JN121115.2','JN121213.2'},{'MP121207.2','MP130705.2','MP130705.4'}}...
    %Set4
    {{'JN121116.2','JN130514.2'},{'MP121210.2','MP130708.2'}}...
    %Set5
    {{'JN121119.2','JN130515.2'},{'MP121211.2','MP130709.2'}}...
    %Set6
    {{'JN121120.2','JN130516.2'},{'MP121212.2','MP130710.3'}}...
    %Set7loo
    {{'JN121121.2','JN130517.2'},{'MP121213.2','MP130711.2'}}...
    };

for setloop=2:size(fList,1);
   for subjloop=1:size(fList{setloop},2) 
       if size(fList{setloop,1}{subjloop},2)>1
           f=1;
           t1=[(2000+str2double(fList{setloop,1}{subjloop}{f}(3:4))) str2double(fList{setloop,1}{subjloop}{f}(5:6)) str2double(fList{setloop,1}{subjloop}{f}(7:8)) 12 0 0];
           f=2;
           t2=[(2000+str2double(fList{setloop,1}{subjloop}{f}(3:4))) str2double(fList{setloop,1}{subjloop}{f}(5:6)) str2double(fList{setloop,1}{subjloop}{f}(7:8)) 12 0 0];
           timeBetween{setloop,1}(subjloop)=((etime(t2,t1)/60)/60)/24;
       end
   end
end

% What Sets Were Shown After Delivering 48 IU of OT? (24 IU/mL) Delivered
% for 10 Minutes, OT=1, SL=0 NaN=No nebulizer
indOT=cell(7,1);
indOT={...
    %Set1 'JN','MP','IW','TT'
    {{nan,0},{nan,0},{nan},{nan}}...
    %Set2 'JN','MP'
    {{1,0},{0,1}}...
    %Set3 'JN','MP'
    {{0,1},{1,0,0}}...
    %Set4 'JN','MP'
    {{1,0},{0,1}}...
    %Set5 'JN','MP'
    {{0,1},{1,0}}...
    %Set6 'JN','MP'
    {{1,0},{0,1}}...
    %Set7 'JN','MP'
    {{0,1},{1,0}}...
    };

save([dirInd 'SSCM_fList_sList_indOT_include',dateSaved,'_',time,'.mat'],'fList','sList','indOT');

% Files With cchgrid calibration
gridCalList={'JN130513.3';'JN130514.2';'JN130515.2';'JN130516.2';'JN130517.2';...
    'MP130704.2';'MP130705.2';'MP130705.4';'MP130708.2';'MP130709.2';'MP130710.3';'MP130711.2'};
% Calibration Files
gridCalFiles{1,1}{1}{2}={'JN130513.1','JN130513.4'};
gridCalFiles{4,1}{1}{2}={'JN130514.1'};
gridCalFiles{5,1}{1}{2}={'JN130515.1','JN130515.3'};
gridCalFiles{6,1}{1}{2}={'JN130516.1','JN130516.3'};
gridCalFiles{7,1}{1}{2}={'JN130517.1','JN130517.3'};

gridCalFiles{1,1}{2}{2}={'MP130703.3','MP130703.5'};
gridCalFiles{2,1}{2}{2}={'MP130704.1','MP130704.3'};
gridCalFiles{3,1}{2}{2}={'MP130705.1','MP130705.5'};
gridCalFiles{3,1}{2}{3}={'MP130705.1','MP130705.5'};
gridCalFiles{4,1}{2}{2}={'MP130708.1','MP130708.3'};
gridCalFiles{5,1}{2}{2}={'MP130709.1','MP130709.3'};
gridCalFiles{6,1}{2}{2}={'MP130710.2','MP130710.4'};
gridCalFiles{7,1}{2}{2}={'MP130711.1','MP130711.3'};



distCutoff = 5;     % If raw calibration data is more than distCutoff DVA away from control point, exclude it from poly transform calculation
sampleTime = 100;  % How many milliseconds to sample before yellow square.
numConditions = 63;  % How our Cortex encodes work right now, this number
                      % can be anything >= your actual condition number, as
                      % long as it is < 1000.  Basically it is 'up to which
                      % condition number do you want to look at?
samprate = 5;   % Samples every samprate ms. of the eyedata, as determined in CORTEX.
calPoints = 63;  % How many / which points to use for the transformation.  How many conditions, starting from 1?

% Combine grid calibration from before & after SSCM?
combineGrid=1;

% Conversion Factor for EOG->DVA
eog2dva=0.008;


% What .itm and .cnd files?
itmFileList={...
    %Set1
    'SSCM90.itm';...
    %Set2
    'SSCMS2.itm';...
    %Set3
    'SSCMS3.itm';...
    %Set4
    'SSCMS4.itm';...
    %Set5
    'SSCMS5.itm';...
    %Set6
    'SSCMS6.itm';...
    %Set7
    'SSCMS7.itm';...
    };
cndFileList={...
    %Set1
    'SSCM90_ERROR.cnd';...
    %Set2
    'SSCMS2.cnd';...
    %Set3
    'SSCMS3.cnd';...
    %Set4
    'SSCMS4.cnd';...
    %Set5
    'SSCMS5.cnd';...
    %Set6
    'SSCMS6.cnd';...
    %Set7
    'SSCMS7.cnd';...
    };
% Initialize array for eye data in dva
eyedat=cell(size(fList,1),1);

% Initialize array for eye data converted to pixel coordinates
eyedatpixlr=cell(size(fList,1),1);

% Initialize array of total Looking Time during scene
ltList=cell(size(fList,1),1);

% Initialize x-y calibration array
xycal=cell(size(fList,1),1);

% Initiliaze Condition List array
cndList=cell(size(fList,1),1);

% Initialize time spent outside image array
timeOutside=cell(size(fList,1),1);

% Initialize Scaling array
xyScale=cell(size(fList,1),1);

%Calibration Points: 9 points in crosshair
% EOG units
calX=[0,0,-375,375,0,0,-750,750,0];
calY=[0,375,0,0,-375,750,0,0,-750];
% figure;
% plot(calX,calY,'.')

%%


fileGrid=1;
tic
% Loop through every Set
for setloop=2:size(fList,1);
%     saveDir=[dirEye num2str(setloop) '\']; 
    
    % Loop through every subject with data for the current Set
    for subjloop=1:size(fList{setloop,1},2)
        % Loop through every Data file in the current Set
        for fileloop=1:size(fList{setloop,1}{subjloop},2);
            
            fileName=fList{setloop,1}{subjloop}{fileloop};
            
            % Select data file
            datfil=[dirDat fileName];
            
            % Select the cnd File for the Set
            if strcmp('JN121108.2',fileName);
                cndFile=[dirCnd,'SSCMS2_ERROR.cnd'];
            elseif strcmp('JN130513.3',fileName);
                cndFile=[dirCnd,'SSCM90.cnd'];
            else
                cndFile=[dirCnd cndFileList{setloop}];
            end
            
            
            %--------------- Select the itm File for the Set-------------------------
            itmFile=[dirItm itmFileList{setloop}];
            
            
            %---------- load itm and cnd files and analyze
            %----------------------------
            itmfil=[];
            [fid,message]=fopen(itmFile, 'r');
            if fid<0
                disp(message);
            else
                while 1
                    tline = fgetl(fid);
                    if ~isempty(itmfil)
                        if length(tline)>size(itmfil,2)
                            tline=tline(1:size(itmfil,2));
                        end
                    end
                    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
                    if ischar(tline)
                        itmfil=[itmfil; tline];
                    else
                        break
                    end
                end
            end
            fclose(fid);
            
            cndfil=[];
            [fid,message]=fopen(cndFile, 'r');
            if fid<0
                disp(message);
            else
                while 1
                    tline = fgetl(fid);
                    if ~isempty(cndfil)
                        if length(tline)>size(cndfil,2)
                            tline=tline(1:size(cndfil,2));
                        end
                    end
                    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
                    if ischar(tline)
                        cndfil=[cndfil; tline];
                    else
                        break
                    end
                end
            end
            fclose(fid);
            
            test0start=strfind(cndfil(1,:),'TEST0');
            test0str=cndfil(:,test0start:test0start+4);
            test0arr=[];
            for k=1:size(test0str,1)
                test0arr=[test0arr; str2double(test0str(k,:))];
            end
            handles.test0arr=test0arr;
            
            cndstart=strfind(cndfil(1,:),'COND#');
            cndstr=cndfil(:,cndstart:cndstart+4);
            cndarr=[];
            for k=1:size(cndstr,1)
                cndarr=[cndarr; str2double(cndstr(k,:))];
            end
            handles.cndarr=cndarr;
            
            trltypstart=strfind(cndfil(1,:),'TRIAL_TYPE');
            trltypstr=cndfil(:,trltypstart:trltypstart+9);
            trltyparr=[];
            for k=1:size(trltypstr,1)
                trltyparr=[trltyparr; str2double(trltypstr(k,:))];
            end
            handles.trltyparr=trltyparr;
            
            itmstart=strfind(itmfil(1,:),'ITEM');
            itmstr=itmfil(:,itmstart:itmstart+3);
            itmarr=[];
            for k=1:size(itmstr,1)
                itmarr=[itmarr; str2double(itmstr(k,:))];
            end
            handles.itmarr=itmarr;
            
            filenamestart=strfind(itmfil(1,:),'------FILENAME------');
            
%             % Repair & Merge MP130705.2 & MP130705.4
%             netdir='S:\Cortex Data\Peepers\';
%             [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata([netdir 'MP130705.2']);
%             
            if ismember(fileName,gridCalList)
                x=cell(1,numConditions);
                y=cell(1,numConditions);
                if combineGrid
                    gridEnd=size(gridCalFiles{setloop,1}{subjloop}{fileloop},2);
                else
                    gridEnd=1;
                end
                hasSamp=cell(1,gridEnd);
                for calFile=1:gridEnd;
                    %Get Data
                    [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata([dirDat gridCalFiles{setloop,1}{subjloop}{fileloop}{calFile}]);
                    
                    event_arr = [event_arr(1:4, :); event_arr(6:end, :)];  % Remove trial-number row of encodes from event_arr because it screws with some later calculations.
                    time_arr = [time_arr(1:4, :); time_arr(6:end, :)];  % Same for the time_arr.
                    calItm=[dirItm 'cchgrid.itm'];
                    stimLocs = readITM(calItm);
                    stimX = stimLocs(1:calPoints, 1);
                    stimY = stimLocs(1:calPoints, 2);
                    
                    
                    calbeg=500;
                    calend=900;
                    % Get eyedata for calibration with clrchng trials
                    numrpt = size(event_arr,2);
                    valrptcnt = 0;
                    clear per clrchgind
                    for rptlop = 1:numrpt
                        if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) <= 1000 + numConditions)) ~=0
                            if size(find(event_arr(:,rptlop) == 200)) ~=0
                                perbegind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                                perendind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                                % Sample calbeg ms after appearance of gray fixspot
                                % until calend ms after appearance of gray fixspot
                                begtimdum = time_arr(perbegind,rptlop)+calbeg;
                                endtimdum = time_arr(perendind,rptlop)+calend;
                                if endtimdum > begtimdum
                                    valrptcnt = valrptcnt + 1;
                                    clrchgind(valrptcnt)=rptlop;
                                    per(valrptcnt).begsmpind = begtimdum;
                                    per(valrptcnt).endsmpind = endtimdum;
                                    per(valrptcnt).begpos = 1;
                                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                    per(valrptcnt).allval = event_arr(:,rptlop);
                                    per(valrptcnt).alltim = time_arr(:,rptlop);
                                end
                            end
                        end
                    end
                    
                    clear cnd cndlst
                    numrpt = size(per,2);
                    for rptlop = 1:numrpt
                        cnd(rptlop)=per(rptlop).cnd;
                        cndList{setloop,1}{1,fileloop}(1,rptlop)=per(rptlop).cnd;
                    end
                    
                    evnnmb=2:2:size(eog_arr,1);
                    oddnmb=1:2:size(eog_arr,1);
                    
                    % Create structures x and y of the corresponding average eye data for each trial
                    % instance (p) of each condition (k)
                    cndlst=unique(cnd);
                    for k=1:length(cndlst)
                        cndind=find(cnd==cndlst(k));
                        allind=clrchgind(cndind);
                        start=size(x{cndlst(k)-1000},2);
                        for p=1:length(allind)
                            x{cndlst(k)-1000}(p+start)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,oddnmb),allind(p))).*eog2dva;
                            y{cndlst(k)-1000}(p+start)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,evnnmb),allind(p))).*eog2dva;
                        end
                    end
                end
                xO=x;
                yO=y;
                
                % Values that are outDVA DVA away from median are
                % removed
                outDVA=1;
                
                for k=1:length(x)
                    if size(x{k},2)>1 && (std(x{k})>outDVA || std(y{k})>outDVA)
                        xMedian=median(x{k});
                        yMedian=median(y{k});
                        xO{k}=x{k}(find(x{k}< xMedian + outDVA & x{k}> xMedian - outDVA ));
                        yO{k}=y{k}(find(x{k}< xMedian + outDVA & x{k}> xMedian - outDVA ));
                        xO{k}=x{k}(find(y{k}< yMedian + outDVA & y{k}> yMedian - outDVA ));
                        yO{k}=y{k}(find(y{k}< yMedian + outDVA & y{k}> yMedian - outDVA ));
                    end
                end
                
                % Create arrays of the average x and y data for each condition across all
                % trial instances of that condition.
                clear meanx meany
                hasSampC=1;
                clear hasSamp
                for k=1:length(xO)
                    if ~isempty(xO{k})
                        meanx(hasSampC)=mean(xO{k});
                        
                        meany(hasSampC)=mean(yO{k});
                        
                        hasSamp(hasSampC)=k;
                        hasSampC=hasSampC+1;
                    end
                end
                control = [stimX(hasSamp,1) stimY(hasSamp,1)];% actual, desired points
                % Begin polynomial transformation
                clear input
                for i = 1:size(meanx,2)
                    input(i, 1) = meanx(i);
                    input(i, 2) = meany(i);
                end
                
                % Compute the polynomial transformation function
                tform = cp2tform(control, input,'polynomial',4);
                tform.forward_fcn = tform.inverse_fcn;
                %             figure
                %             for j = 1:size(x, 2)
                %                     scatter(x{j}, y{j}, 7, 'r', 'filled');
                %                     hold all
                %             end
                %
                %             % Create figure, plot the raw data with the forward function applied
                %             figure
                %             for j = 1:size(x, 2)
                %                     [x2{j} y2{j}] = tformfwd(tform, xO{j}, yO{j});
                %                     scatter(x2{j}, y2{j}, 7, 'r', 'filled');
                %                     hold all
                %             end
                % %             % Apply the forward function to the mean eye data and add them to the figure.
                %             [meanx2 meany2] = tformfwd(tform, meanx, meany);
                %             scatter(meanx2, meany2, 'g', 'filled');
                %             axis([-20 20 -20 20])
                %             set(gca, 'Color', 'k');
                %
                %             %  Add the actual stimuli locations to the plot
                %             scatter(stimX, stimY, 'y', 'filled', 'marker', 's');
                %             title(['Setloop=',num2str(setloop),' File=',fileName]);
            else
                clear x y
%                 clear time_arr event_arr eog_arr header trialcount
                % Get Data
                [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datfil);
                
                calbeg=500;
                calend=900;
                % Get eyedata for calibration with clrchng trials
                numrpt = size(event_arr,2);
                valrptcnt = 0;
                clear per clrchgind
                for rptlop = 1:numrpt
                    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 1010)) ~=0
                        if size(find(event_arr(:,rptlop) == 200)) ~=0
                            perbegind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                            perendind = find(event_arr(:,rptlop) == 23);%color change spot appears in gray
                            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                            % Sample calbeg ms after appearance of gray fixspot
                            % until calend ms after appearance of gray fixspot
                            begtimdum = time_arr(perbegind,rptlop)+calbeg;
                            endtimdum = time_arr(perendind,rptlop)+calend;
                            if endtimdum > begtimdum
                                valrptcnt = valrptcnt + 1;
                                clrchgind(valrptcnt)=rptlop;
                                per(valrptcnt).begsmpind = begtimdum;
                                per(valrptcnt).endsmpind = endtimdum;
                                per(valrptcnt).begpos = 1;
                                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                                per(valrptcnt).allval = event_arr(:,rptlop);
                                per(valrptcnt).alltim = time_arr(:,rptlop);
                            end
                        end
                    end
                end
                
                clear cnd
                numrpt = size(per,2);
                for rptlop = 1:numrpt
                    cnd(rptlop)=per(rptlop).cnd;
%                     cndList{setloop,1}{1,fileloop}(1,rptlop)=per(rptlop).cnd;
                end
                
                evnnmb=2:2:size(eog_arr,1);
                oddnmb=1:2:size(eog_arr,1);
                
                clear x y
                cndlst=unique(cnd);
                for k=1:length(cndlst)
                    cndind=find(cnd==cndlst(k));
                    allind=clrchgind(cndind);
                    for p=1:length(allind)
                        x{k}(p)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,oddnmb),allind(p))).*eog2dva;
                        y{k}(p)=mean(eog_arr(intersect(floor(((per(cndind(p)).begsmpind-1000)/samprate)*2):(floor((per(cndind(p)).endsmpind-1000)/samprate))*2,evnnmb),allind(p))).*eog2dva;
                    end
                end
                
                if strmatch('IW120425.1',fileName);
                    %Correct for error in cnd 2 and rename of cnd 5 to cnd 2 in IW120425.1
                    %Remove Condition 2 & 3
                    x={x{1,1},x{1,4},x{1,5},x{1,6},x{1,7},x{1,8},x{1,9}};
                    y={y{1,1},y{1,4},y{1,5},y{1,6},y{1,7},y{1,8},y{1,9}};
                    %Change Calibration Points: Remove Cnd 2 and put 5 at 0,6 degrees
                    calX=[0,  3,  0,   0,  -6, 6,  0];
                    calY=[0,  0,   6, 6,   0,   0,  -6];
                elseif strmatch('JN121108.2',fileName);
                    %Remove Condition 3
                    x={x{1,1},x{1,2},x{1,4},x{1,5},x{1,6},x{1,7},x{1,8},x{1,9}};
                    y={y{1,1},y{1,2},y{1,4},y{1,5},y{1,6},y{1,7},y{1,8},y{1,9}};
                    %Change Calibration Points: Put Cnd 5 at 0,6 degrees
                    calX=[0,0,3,0,0,-6,6,0];
                    calY=[0,3,0,6,6,0,0,-6];
                elseif strmatch('JN121212.2',fileName);
                    %Remove Condition 3
                    x={x{1,1},x{1,2},x{1,4},x{1,5},x{1,6},x{1,7},x{1,8},x{1,9}};
                    y={y{1,1},y{1,2},y{1,4},y{1,5},y{1,6},y{1,7},y{1,8},y{1,9}};
                    %Change Calibration Points: Put Cnd 5 at 0,6 degrees
                    calX=[0,0,3,0,0,-6,6,0];
                    calY=[0,3,0,6,6,0,0,-6];
                else
                    calX=[0,0,-3,3,0,0,-6,6,0];
                    calY=[0,3,0,0,-3,6,0,0,-6];
                end
                clear xO yO
                xO=x;
                yO=y;
                
                % Values that are outDVA DVA away from median are
                % removed
                outDVA=1;
                
                for k=1:length(x)
                    if size(x{k},2)>1 && (std(x{k})>outDVA || std(y{k})>outDVA)
                        xMedian=median(x{k});
                        yMedian=median(y{k});
                        xO{k}=x{k}(find(x{k}< xMedian + outDVA & x{k}> xMedian - outDVA ));
                        yO{k}=y{k}(find(x{k}< xMedian + outDVA & x{k}> xMedian - outDVA ));
                        xO{k}=x{k}(find(y{k}< yMedian + outDVA & y{k}> yMedian - outDVA ));
                        yO{k}=y{k}(find(y{k}< yMedian + outDVA & y{k}> yMedian - outDVA ));
                    end
                end
                
                % Create arrays of the average x and y data for each condition across all
                % trial instances of that condition.
                clear meanx meany
                for k=1:length(xO)
                    meanx(k)=mean(xO{k});
                    meany(k)=mean(yO{k});
                end
                
                clear input
                for i = 1:size(meanx,2)
                    input(i, 1) = meanx(i);
                    input(i, 2) = meany(i);
                end
                
                control = [calX' calY'];% actual, desired points
                tform = cp2tform(control, input,'affine');
                tform.forward_fcn = tform.inverse_fcn;
                %             % Create figure, plot the raw data with the forward function applied
                %             figure
                %             for j = 1:size(x, 2)
                %                 for i = 1:size(x{j},2)
                %                     [x2{j}( :, i) y2{j}( :, i)] = tformfwd(tform, xO{j}( :, i), yO{j}( :, i));
                %                     scatter(x2{j}( : , i), y2{j}( :, i), 7, 'r', 'filled');
                %                     hold all
                %                 end
                %             end
                %             % Apply the forward function to the mean eye data and add them to the figure.
                %             [meanx2 meany2] = tformfwd(tform, meanx, meany);
                %             scatter(meanx2, meany2, 'g', 'filled');
                %             axis([-20 20 -20 20])
                %             set(gca, 'Color', 'k');
                %
                %             %  Add the actual stimuli locations to the plot
                %             scatter(calX, calY, 'y', 'filled', 'marker', 's');
                %             title(['Setloop=',num2str(setloop),' File=',fileName]);
                
            end
            if ismember(fileName,gridCalList)
                clear time_arr event_arr eog_arr header trialcount
                % Get SSCM Data
                [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datfil);
            end
            
            % Get eye data for stimulus presentation
            numrpt = size(event_arr,2);
            valrptcnt = 0;
            clear per vpcind
            new_eog_arr=[];
            for rptlop = 1:numrpt
                % Find all stimulus presentation trials (cnd>=1010)
                if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
                    if size(find(event_arr(:,rptlop) == 200)) ~=0
                        % 23 instead of 24 because 23 is stim onset
                        perbegind = find(event_arr(:,rptlop) == 23,1,'first');
                        perendind = find(event_arr(:,rptlop) == 24,1,'first');
                        cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                        blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                        begtimdum = time_arr(perbegind,rptlop);
                        endtimdum = time_arr(perendind,rptlop);
                        if endtimdum > begtimdum
                            valrptcnt = valrptcnt + 1;
                            vpcind(valrptcnt)=rptlop;
                            per(valrptcnt).begsmpind = begtimdum;
                            per(valrptcnt).endsmpind = endtimdum;
                            per(valrptcnt).begpos = 1;
                            per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                            per(valrptcnt).blk = event_arr(blknumind,rptlop);
                            per(valrptcnt).allval = event_arr(:,rptlop);
                            per(valrptcnt).alltim = time_arr(:,rptlop);
                            %                         new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                            new_eog_arr=[new_eog_arr,eog_arr(:,rptlop)];
                        end
                    end
                end
            end
            
            clear cnd cndlst
            numrpt = size(per,2);
            for rptlop = 1:numrpt
                cnd(rptlop)=per(rptlop).cnd;
%                 cndList{setloop,1}{1,fileloop}(1,rptlop)=per(rptlop).cnd;
                startInd=find(expectedCnd==cnd(1));
            end
            
            % Check that all the expected conditions were run
            cndShown(1,1:size(expectedCnd,2))=false;
            for k=1:size(cnd,2)
                cndShown(k)=cnd(k)==expectedCnd(k);
            end
            
            for trlop=1:size(per,2)
%                 if cndShown(trlop);
                    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
                    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
                    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
                    picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to iscan start
                    picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to iscan start
                    
                    if picend/5<=length(horeog)
                        xDat=(horeog(round(picstart/5):floor(picend/5))).*eog2dva';
                        yDat=(vrteog(round(picstart/5):floor(picend/5))).*eog2dva';
                        [eyedat{setloop,1}{subjloop}{fileloop}{trlop,1}(1,:),eyedat{setloop,1}{subjloop}{fileloop}{trlop,1}(2,:)] = tformfwd(tform,xDat,yDat);
                    else
                        eyedat{setloop,1}{subjloop}{fileloop}{trlop,1}(1,:)=nan;
                        eyedat{setloop,1}{subjloop}{fileloop}{trlop,1}(2,:)=nan;
                    end
                
%                 end
            end
            
            if strcmp(fileName,'MP130705.2')
                eyedat{setloop,1}{subjloop}{fileloop}=eyedat{setloop,1}{subjloop}{fileloop}(cndShown,1);
            end
            
            %---remove values outside of 800x600 pixel screen---%
            badxCell=cell(size(eyedat{setloop,1}{subjloop}{fileloop},1),1);
            badyCell=cell(size(eyedat{setloop,1}{subjloop}{fileloop},1),1);
            for i = 1:size(eyedat{setloop,1}{subjloop}{fileloop},1);
                x = (eyedat{setloop,1}{subjloop}{fileloop}{i,1}(1,:)*24)+400;
                y = (eyedat{setloop,1}{subjloop}{fileloop}{i,1}(2,:)*24)+300;
                badxCell{i} = find(x < 1 | x > 800);
                badx = find(x < 1 | x > 800);
                x(badx) = []; y(badx) = [];
                badyCell{i} = find(y < 1 | y > 600);
                bady = find(y < 1 | y > 600);
                x(bady) = []; y(bady) = [];
                x = (x-400)/24; y = (y-300)/24;
                eyedat{setloop,1}{subjloop}{fileloop}{i,1} = [x;y];
            end
            timeOutside{setloop,1}{subjloop}{fileloop}=zeros(size(badxCell));
            for k=1:length(badxCell);
                timeOutside{setloop,1}{subjloop}{fileloop}(k,1)=length(badxCell{k,1})*5;
            end
            
            clear lt cnd
            numrpt = size(per,2);
            for rptlop = 1:numrpt
                ltList{setloop,1}{subjloop}{fileloop}{rptlop,1}=per(rptlop).endsmpind-per(rptlop).begsmpind;
            end
            
            %Convert Eye data units from degrees of visual angle to pixel coordinates
            for k=1:size(eyedat{setloop,1}{subjloop}{fileloop},1)
                eyedatpixlr{setloop,1}{subjloop}{fileloop}{k,1}(1,:)=round(((eyedat{setloop,1}{subjloop}{fileloop}{k,1}(1,:)*24)+400));
                eyedatpixlr{setloop,1}{subjloop}{fileloop}{k,1}(2,:)=round(((eyedat{setloop,1}{subjloop}{fileloop}{k,1}(2,:)*24)+300));
            end
            
            % MP stopped during 'MP130705.2' and JAS hit escape too many times and had to start a new file, and task
            % was mistakenly resumed in wrong place and scenes 62, 63, and 64 were skipped
            if strcmp(fileName,'MP130705.4')
                % Last condition shown was the first scene (task restarted
                % from beginning because the number of conditions shown
                % wasn't changed before starting 'MP130705.4' 
                eyedat{setloop,1}{subjloop}{2}(startInd:180)=eyedat{setloop,1}{subjloop}{3}(1:end-1);
                eyedatpixlr{setloop,1}{subjloop}{2}(startInd:180)=eyedatpixlr{setloop,1}{subjloop}{3}(1:end-1);
                timeOutside{setloop,1}{subjloop}{2}(startInd:180)=timeOutside{setloop,1}{subjloop}{3}(1:end-1);
                timeOutside{setloop,1}{subjloop}{2}(123:128,1)=nan;
                
                % Get rid of 'MP130705.4' data
                eyedat{setloop,1}{subjloop}=eyedat{setloop,1}{subjloop}(1:2);
                eyedatpixlr{setloop,1}{subjloop}=eyedatpixlr{setloop,1}{subjloop}(1:2);
                timeOutside{setloop,1}{subjloop}=timeOutside{setloop,1}{subjloop}(1:2);
            end
        end
    end
%     if ~exist(saveDir,'dir');
%         mkdir(saveDir);
%     end
%     save([saveDir 'SSCM_eyedat_S',num2str(setloop),'_',dateSaved,'_',time,'.mat'],'indOT','fList','sList','eyedat','itmFileList','cndFileList');
%     clear eyedat
%     save([saveDir 'SSCM_eyedatpixlr_S',num2str(setloop),'_',dateSaved,'_',time,'.mat'],'indOT','fList','sList','eyedatpixlr','itmFileList','cndFileList');
%     clear eyedatpixlr
end
toc

%%
dirVar='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';

dateSaved=date;
time=clock;
time=[num2str(time(4)),'h',num2str(time(5)),'m'];
% Save eyedatpixlr only
save([vardir 'SSCM_EyeDataPixl_',dateSaved,'_',time,'.mat'],'indOT','fList','sList','eyedatpixlr','itmFileList','cndFileList');
% Save eyedat only
save([vardir 'SSCM_eyedat_',dateSaved,'_',time,'.mat'],'indOT','fList','sList','eyedat','itmFileList','cndFileList');
% Save timeOutside only
save([vardir 'SSCM_timeOutside_',dateSaved,'_',time,'.mat'],'indOT','fList','sList','timeOutside','itmFileList','cndFileList');

%% Having problems getting data from data files
vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
eyeDir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\eyedat\';
% Load Eye Position During Stimulus Presentation for Sets Specified Above
% Recalibration with cp2tform affine transformation for files with 9-point
% cross and polynomial transformation for files with 63-point grid 

for setloop=1:7
    load([vardir 'SSCM_eyedat_12-Jul-2013_17h11m.mat']);
    for k=1:size(eyedat,1)
        if k~=setloop
            eyedat{k,1}=[];
        end
    end
    saveDir=[eyeDir num2str(setloop) '\'];
    if ~exist(saveDir,'dir');
        mkdir(saveDir);
    end
    save([saveDir 'SSCM_eyedat_',num2str(setloop),'.mat'],'indOT','fList','sList','eyedat','itmFileList','cndFileList');
end


%% Plot Time Outside
% Plot Each Session Separately
clear evnnmb oddnmb
evnnmb=2:2:T;
oddnmb=1:2:T;


z=1;
ot=[1,1];
sl=[1,1];
clear timeOutsideAll D
for setloop=2:length(timeOutside);
    for subjloop=1:length(timeOutside{setloop,1})
        for fileloop=1:length(timeOutside{setloop,1}{subjloop})
%             figure;
%             plot(timeOutside{setloop,1}{subjloop}{fileloop});
%             title(fList{setloop,1}{subjloop}{fileloop});
            z=z+1;
            timeOutsideAll(z,:)=timeOutside{setloop,1}{subjloop}{fileloop};
            
            if indOT{setloop,1}{subjloop}{fileloop};
                D{2,fileloop}{1}(ot(fileloop),:)=timeOutside{setloop,1}{subjloop}{fileloop}(oddnmb);
                D{2,fileloop}{2}(ot(fileloop),:)=timeOutside{setloop,1}{subjloop}{fileloop}(evnnmb);
                ot(fileloop)=ot(fileloop)+1;
            else
                D{1,fileloop}{1}(sl(fileloop),:)=timeOutside{setloop,1}{subjloop}{fileloop}(oddnmb);
                D{1,fileloop}{2}(sl(fileloop),:)=timeOutside{setloop,1}{subjloop}{fileloop}(evnnmb);
                sl(fileloop)=sl(fileloop)+1;
            end
        end
    end
end

clear group anovaD
c=1;
for drug=1:2
    for pres=1:2
        for trl=1:2
            for fileloop=1:size(D{pres,drug}{trl},1)
                dat=D{pres,drug}{trl}(fileloop,:)>0;
                anovaD(c:c+(length(D{pres,drug}{trl})-1))=D{pres,drug}{trl}(fileloop,:);
                group{1,1}(c:c+(length(D{pres,drug}{trl}(fileloop,:))-1))=drug;
                group{2,1}(c:c+(length(D{pres,drug}{trl}(fileloop,:))-1))=pres;
                group{3,1}(c:c+(length(D{pres,drug}{trl}(fileloop,:))-1))=trl;
                c=c+length(D{pres,drug}{trl}(fileloop,:));
            end
        end
    end
end
[p,table,stats]=anovan(anovaD,group,'varnames',{'Drug','Presentation','Trial'},'model','interaction','display','on')


[P,~,stats] = anovan(aDat,{gGaze,gDrug},'varnames',{'Gaze','Drug'},'model','interaction','display','on');
[compMat,meanMat]=multcompare(stats,'display','on','dimension',[1 2 3]);
mean(anovaD(group{3,1}==2))

[h,p]=ttest2(sum(D{2,1},2),sum(D{2,2},2));

figure;
plot([(mean(sum(D{1,1},2))/1000),(mean(sum(D{2,1},2))/1000)],'b')
hold on;
plot([(mean(sum(D{1,2},2))/1000),(mean(sum(D{2,2},2))/1000)],'r')
ylim([0 15])
channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% Plot Each Session Separately
figure;
z=1;
for s=1:length(timeOutside);
    for f=1:length(timeOutside{s,1})
        plot(timeOutside{s,1}{1,f});
        hold on;
        z=z+1;
        timeOutsideAll(z,:)=timeOutside{s,1}{1,f};
    end
end
%Plot Average
smval=10;
fsamp=0;
figure;
dofill(1:180,timeOutsideAll,'b',1,smval,0,0,1,fsamp,0)
figure;
plot(nanmean(timeOutsideAll));