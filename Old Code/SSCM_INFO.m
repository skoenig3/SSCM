%% SSCM INFO

% Information about filenames, drug treatment, and other SSCM methods

%% Location of Saved Data/Indices
vardir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/';
inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';


%% TRIAL SETTINGS/CONVERSION FACTORS

% Conversion Factor for EOG->DVA
eog2dva=0.008;

% Number of Trials?
T=180;

%ISCAN Sampling Frequency
sampRate=200;

%% SUBJECTS

% What Monkeys?
MonkeyList={...
    %Set1
    {'JN','MP','IW','TT','JN'};...
    %Set2
    {'JN','MP','JN'};...
    %Set3
    {'JN','MP','JN'}...
    %Set4
    {'JN','MP','JN'}...
    %Set5
    {'JN','MP','JN'}...
    %Set6
    {'JN','MP','JN'}...
    %Set7
    {'JN','MP','JN'}...
    };

%% FILENAMES

% What data files?
Sets={...
    %Set1
    {'JN120717.4','MP120716.1','IW120425.1','TT120717.5','JN130513.3'};...
    %Set2
    {'JN121108.2','MP121206.2','JN121212.2'};...
    %Set3
    {'JN121115.2','MP121207.2','JN121213.2'}...
    %Set4
    {'JN121116.2','MP121210.2','JN130514.2'}...
    %Set5
    {'JN121119.2','MP121211.2','JN130515.2'}...
    %Set6
    {'JN121120.2','MP121212.2','JN130516.2'}...
    %Set7
    {'JN121121.2','MP121213.2','JN130517.2'}...
    };

presNum={...
    %Set 1
    {1,1,1,1,2};...
    %Set2
    {1,1,2};...
    %Set3
    {1,1,2}...
    %Set4
    {1,1,2}...
    %Set5
    {1,1,2}...
    %Set6
    {1,1,2}...
    %Set7
    {1,1,2}...
    };

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

% Files separated by set and subject
fList={...
    %Set1
    {{'JN120717.4','JN130513.3'},{'MP120716.1'},{'IW120425.1'},{'TT120717.5'}}...
    %Set2
    {{'JN121108.2','JN121212.2'},{'MP121206.2'}}...
    %Set3
    {{'JN121115.2','JN121213.2'},{'MP121207.2'}}...
    %Set4
    {{'JN121116.2','JN130514.2'},{'MP121210.2'}}...
    %Set5
    {{'JN121119.2','JN130515.2'},{'MP121211.2'}}...
    %Set6
    {{'JN121120.2','JN130516.2'},{'MP121212.2'}}...
    %Set7
    {{'JN121121.2','JN130517.2'},{'MP121213.2'}}...
    };


%% INDICES

% % Find PRESENTATION # Automatically
% for setloop=1:size(Sets,1)
%     s=1;
%     for fileloop=1:size(Sets{setloop,1},2)
%         if fileloop==1
%             sList{s}=Sets{setloop,1}{fileloop}(1:2);
%             sListC(s)=1;
%         else
%             for list=1:size(sList,2);
%                 if strcmp(Sets{setloop,1}{fileloop}(1:2),sList{list})
%                     
%                 end
%             end
%         end
%     end
%     presNum{setloop,1}(fileloop)=1;
%     
% end





% What Sets Were Shown After Delivering 48 IU of OT? (24 IU/mL) Delivered
% for 10 Minutes, OT=1, SL=0
indOT=cell(7,1);
indOT{1,1}=logical([0,0,0,0,0]);
%Set 2 JN1,MP, JN2
indOT{2,1}=logical([1,0,0]);
%Set 3 JN, MP, JN2
indOT{3,1}=logical([0,1,1]);
%Set 4 JN, MP
indOT{4,1}=logical([1,0,0]);
%Set 5 JN, MP
indOT{5,1}=logical([0,1,1]);
%Set 6 JN, MP
indOT{6,1}=logical([1,0,0]);
%Set 7 JN, MP
indOT{7,1}=logical([0,1,1]);


%% CORTEX FILES

% Files With cchgrid calibration
gridCalList={'JN130513.3';'JN130514.2';'JN130515.2';'JN130516.2';'JN130517.2'};
% Calibration Files
gridCalFiles{1,1}{1,5}={'JN130513.1','JN130513.4'};
gridCalFiles{4,1}{1,3}={'JN130514.1'};
gridCalFiles{5,1}{1,3}={'JN130515.1','JN130515.3'};
gridCalFiles{6,1}{1,3}={'JN130516.1','JN130516.3'};
gridCalFiles{7,1}{1,3}={'JN130517.1','JN130517.3'};

% What .itm and .cnd files?
datdir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Data Files/';
itmdir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Itm Files/';
cnddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Cnd Files/';
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