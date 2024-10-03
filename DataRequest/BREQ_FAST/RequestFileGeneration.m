%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate request data files via email (BREQ_FAST)
%% The formate of the request file can be found in http://ds.iris.edu/ds/nodes/dmc/manuals/breq_fast/
%%
%% Note: Station and event Catalog shoud be prepared, the format is like below 
% Station Catalog:
%   Network Station    Lat.      Lon.   Start_time   End_time   Channel1 Channel2 Channel3
%    CN     A11     47.24310 70.196800 1994/10/16   2017/10/16   HHE     HHN      HHZ
%     ...
%  Event Catalog:
%  Date         Time        Lat.     Lon.     Depth  Magnitude_tag  Magnitude
%  2017-10-24 16:12:14.880 15.2536 -94.0342   44.91      mb            5.0
%     ...
%
%%    Cong Li     04/15/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
Porjdir=['/Users/congli/Project/ReceiverFunc_Appalachians/'];
Catadir=['Figure/ShareData/catalog/'];
Workdir=['DataRequest/EQs_request/'];
[nws, stas, slat, slon, sdate, edate, schanel1,schanel2,schanel3] = textread([Porjdir Catadir 'StationCatalog.txt'],'%s %s %f %f %s %s %s %s %s');
[date,time,lat,lon,dep,magtag,mag] = textread([Porjdir Catadir 'EQs_MB5.0_2000_2016.events.txt'],'%s %s %f %f %f %s %f');
nevents = length(date);
for nfiles = 1:length(slon)% Change the format of station and network name to characteristic type    
sta = char(stas(nfiles));
nw = char(nws(nfiles));    
fname=[Porjdir Workdir 'Events_' char(sta) '_' num2str(nfiles)];    
fid = fopen(fname,'w');

fprintf(fid, '%s \n', '.NAME CongLI ');
fprintf(fid, '%s \n', '.INST IRIS DMC ');
fprintf(fid, '%s \n', '.EMAIL leecong87@gmail.com ');
fprintf(fid, '%s \n', ['.LABEL ' fname  ]);

fprintf(fid, '%s \n', '.PHONE 413-577-1250 ');
fprintf(fid, '%s \n', '.MEDIA FTP ');
fprintf(fid, '%s \n', '.ALTERNATE MEDIA ');
fprintf(fid, '%s \n', '.ALTERNATE MEDIA ');
fprintf(fid, '%s \n', '.FROM_JWEED ');
fprintf(fid, '%s \n', '.END ');
idx1 = find(datenum(date) >= datenum(sdate{nfiles})& datenum(date) <= datenum(edate{nfiles}));   
for nn = idx1(1):idx1(length(idx1))
    
%  year = date{nn}(1:4);
%  mon = date{nn}(6:7);
%  day = date{nn}(9:10);
%  hr = time{nn}(1:2);
%  min = time{nn}(4:5);
%  sec = time{nn}(7:8);
%  fprintf(fid, ' %s %s %d %d %d %d %d %d %d %d %d %d %d %d %d %s %s %s %s \n', sta,nw,str2num(year),str2num(mon),str2num(day),str2num(hr)-1,str2num(min), str2num(sec), str2num(year),str2num(mon),str2num(day),str2num(hr)+1,str2num(min), str2num(sec), 3, 'BHE','BHN','BHZ', '*' );

year = str2num(date{nn}(1:4));
mon = str2num(date{nn}(6:7));
day = str2num(date{nn}(9:10));
hr = str2num(time{nn}(1:2));
min = str2num(time{nn}(4:5));
sec = str2num(time{nn}(7:8));

hr0=hr-1; hr1=hr+1; % Extract the data for 2 hours( 1 hour before the earthquake happen, and 1 hour after the earthquake happen)
day0=day; day1=day;
year0=year;
year1=year;
mon0=mon;
mon1=mon;
if hr==0         %  Avoid the wrong date
    day0=day-1;
    hr0=23;
end 
if hr==23        %  Avoid the wrong date
    day1=day+1;
    hr1=0;
end
if day0==0   
      if mon==2||mon==4||mon==6||mon==8||mon==9||mon==11
          mon0=mon-1;
          day0=31;
      elseif mon==5||mon==7||mon==10||mon==12
          mon0=mon-1;
          day0=30;
      elseif mon==1
          mon0=12;
          day0=31;
          year0=year-1;
      elseif mon==3
          if mod(year,4)==0
              mon0=mon-1;
              day0=29;
          else 
              mon0=mon-1;
              day0=28;
          end
      end
end
if mon==1||mon==3||mon==5||mon==7||mon==8||mon==10
    if day1>31
          mon1=mon+1;
          day1=1;
    end
elseif mon==4||mon==6||mon==9||mon==11
    if day1>30
          mon1=mon+1;
          day1=1;
    end
elseif mon==12
    if day1>31
          mon1=1;
          day1=1;
          year1=year+1;
    end
elseif mon==2
    if mod(year,4)==0
        if day1>29
              mon1=mon+1;
              day1=1;
        end
    else 
        if day1>28
              mon1=mon+1;
              day1=1;
        end
    end
end
if hr0<0 | hr1>=24
    continue
end
if day0<0 | day1>31
    continue
end
        
fprintf(fid, ' %s %s %d %d %d %d %d %d %d %d %d %d %d %d %d %s %s %s %s \n', sta,nw,year0,mon0,day0,hr0,min,sec,year1,mon1,day1,hr1,min,sec, 3, schanel1{nfiles},schanel2{nfiles},schanel3{nfiles}, '*' );
end
fclose(fid);
end
