itmdir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Itm Files/';
itmFile=[itmdir 'CCH109.itm'];

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

itmstart=strfind(itmfil(1,:),'ITEM');
itmstr=itmfil(:,itmstart:itmstart+3);
itmarr=[];
for k=1:size(itmstr,1)
    itmarr=[itmarr; str2double(itmstr(k,:))];
end
handles.itmarr=itmarr;

xstart=strfind(itmfil(1,:),'CENTERX');
xstr=itmfil(:,xstart:xstart+6);
xarr=[];
for k=1:size(xstr,1)
    xarr=[xarr; str2double(xstr(k,:))];
end
handles.xarr=xarr;

ystart=strfind(itmfil(1,:),'CENTERY');
ystr=itmfil(:,ystart:ystart+6);
yarr=[];
for k=1:size(ystr,1)
    yarr=[yarr; str2double(ystr(k,:))];
end
handles.yarr=yarr;
%% Change Item Numbers

for k=10:99
    newItm=num2str(k);
    itmfil((k+5),itmstart+2:itmstart+(1+size(newItm,2)))=newItm;
end

for k=100:329
    newItm=num2str(k);
    itmfil((k+5),itmstart+2:itmstart+(1+size(newItm,2)))=newItm;
end

%% Change Center X & Center Y
spacing=3;
xRange=-15:spacing:15;
yRange=-12:spacing:12;
for y=1:size(yRange,2)
    for x=1:size(xRange,2)
        if xRange(1,x)<0;
            valStart=2;
        else
            valStart=3;
        end
        xStr=[num2str(xRange(1,x)),'.00'];
        strSize=size(xStr,2);
        itmfil(((x*y)+8),xstart+valStart:(xstart+valStart+strSize-1))=xStr;
    end
end
