%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Extraction RFs in each station from Antelope Database, and save the data in mat format
%    Note: Antelope Database should be prepared    
%   
%% History 
% 11/29/2017
% 18/06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
addpath('/Users/congli/Software/Mat_Lib');
AntelopePath   % set Antelope path
run('/opt/antelope/5.7/setup.m');
samprate=40;                % desired sample rate for data
%% ===== extract data from Antelope database =====
 dbname = {'TraceDeconWif'}; %{'SelectRFs_6000cri'};%
 savepath=['/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data'];
 dir = ['/Users/congli/Project/ReceiverFunc_Appalachians/Database/NewData_Appalachians_RFs_2010_2017/' char(dbname)  ];
 db = dbopen( dir, 'r' );
 dbar = dblookup_table( db, 'arrival' );
 dbas = dblookup_table( db, 'assoc' );
 dbas = dbsubset( dbas, 'delta >= 30');
 dbas = dbsubset( dbas, 'delta <= 95');
 dbwf = dblookup_table( db, 'wfdisc' );
 dbar = dbjoin( dbar, dbas );
 dbor = dblookup_table( db, 'origin');
 dbor = dbsubset( dbor, 'mb >= 5.0');
 dbar = dbjoin( dbar,dbor);
 dbst = dblookup_table( db, 'site');
 dbnw = dblookup_table( db,'snetsta');
 dbchan=dblookup_table(db,'sitechan');
 dbori= dblookup_table( db, 'origin' );
 dbchan = dbjoin(dbchan,dbori);
 dbar = dbjoin(dbar,dbwf);
 [artime,arsta,bazimuth,azimuth,arid,orid,edep,lat_ori,lon_ori,del,Mb] = ...
dbgetv(dbar,'time','sta','seaz','esaz','arid','orid','depth','lat','lon','delta','mb');
%  strtime(artime)
%uarsta=unique(arsta);
uorid = unique(orid);
ind = 1;
display(['number of station ' num2str(length(uorid))]);
preA=-30;
postA=120;
%% =========================================================================
for evt = 1577:length(uorid)
%for evt = 1
   run('/opt/antelope/5.7/setup.m');
   fprintf('Origin Id %d: %d\n',evt, uorid(evt));
      
   slowP=[]; ymdtime={}; istn={};
   iorid=[]; olat=[]; olon=[]; slat=[]; slon=[]; baz=[]; atime=[]; dist=[];
   arrivaltime=[]; dep=[];MagB=[];Ehang=[];Nhang=[];Zhang=[];delta=[];
   RF_Zcomp={};RF_Rcomp={};RF_Tcomp={};save_tick=0;sampleo=[];arrtime=[];
   index11 = find(orid == uorid(evt));  % changed from the origonal program get rid of repeated stations for each event
   [trash, index22] = unique(arsta(index11));
   index = index11(index22);
   nStn = length(index);     % number of stations for this event
   count = 0;
   for ind =1:nStn   %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   iar = index(ind);
   [slowness, phasenames] = arr_slowness( del(iar), edep(iar));
%    slowP(ind) = slowness(1,1);
   numPhase=length(slowness);
     for phaseE=1:numPhase       
           pp=cell2mat(phasenames(phaseE));
         if length(pp)==5 & pp(1:5)=='PKPdf';
           phaseNum=phaseE;
          elseif length(pp)==1 && pp(1)=='P'
          phaseNum=phaseE;    
        end
     end
  slowP=slowness(phaseNum);  % assumes phase is P or PKIKP
   dbsub = dbsubset( dbwf, ['sta == "' arsta{iar} '"'] );
%    dbsub  = dbsort(dbsub1, 'chan','time','dbSORT_UNIQUE');
   tr = trload_css(dbsub, artime(iar)+preA, artime(iar)+postA);% <---
   trsplice( tr, 0.5 );
   trapply_calib(tr);
   nrecs = dbnrecs(tr);

    if ( dbnrecs(tr) ~= 3 )
        disp('Not 3 components ...'); %% if there is a dupulited data, error happened, need to fix it later
         continue
    else
        tr.record = 0;
        I0 = trextract_data(tr);
        tr.record = 1;
        I1 = trextract_data(tr);
        tr.record = 2;
        I2 = trextract_data(tr);
        sps = dbgetv( tr, 'samprate' );
        if sps~=round(sps)
            sps=round(sps);
        end
        titt=0;
        if (length(I0)~=length(I1))
             titt=1; 
            continue;
        end
        if ~isempty(find(isnan(I2))), continue, end
    end
  
   nrecs = dbnrecs(tr);
   dbsub = dbsubset( dbst, ['sta == "' arsta{iar} '"'] );
   [slt,sln,elv]  = dbgetv( dbsub, 'lat','lon','elev');
   
   dbsub = dbsubset( dbnw, ['sta == "' arsta{iar} '"'] );
   nwsta  = dbgetv( dbsub, 'snet');
   
   dbsub = dbsubset( dbchan, ['sta == "' arsta{iar} '" && ' sprintf(' orid ==%f ', uorid(evt)) '']);
   [hang,staname] = dbgetv( dbsub, 'hang','sta');
%    if (hang(1)~=90)
%        ;
%    end
%    if isempty(char(arsta(iar))), continue, end
%    
%    if abs(length(HZ)-(postA-preA)*sps)>=2, continue, end
%    if abs(length(H1)-(postA-preA)*sps)>=2, continue, end
   
   count = count+1;
   iorid = orid(iar);
   istn= {[char(nwsta) '.' char(arsta(iar))]};
   olat  = lat_ori(iar);
   olon  = lon_ori(iar);
   slat  = mean(slt);
   slon  = mean(sln);
   selv  = mean(elv);
 if (size(hang,1) == 3)
     Ehang=hang(1);
     Nhang=hang(2);
     Zhang=hang(3);
 elseif (size(hang,1)==2)
     Ehang=hang(1);
     Nhang=hang(2);
     Zhang=9999;    % 9999 means unreliable value
 elseif (size(hang,1)==1)
     Ehang=hang(1);
     Nhang=9999;
     Zhang=9999;
 else
     Ehang=9999;
     Nhang=9999;
     Zhang=9999;
 end

   baz   = bazimuth(iar);
   delta  = del(iar);
   dep = edep(iar);
   MagB = Mb(iar);
   atime = artime(iar);
   ymdtime = strydtime(atime); % convert arrival time to hhmmss
   aa = strydtime(atime); % convert arrival time to hhmmss
   arrt = aa(18:end);
   arrtime = str2num(arrt(1:2))*3600 + ...
                      str2num(arrt(4:5))*60 + ...
                      str2num(arrt(7:8)) + ...
                      str2num(arrt(10:end))*0.001;
   
   RF_Zcomp = I2(:);
   RF_Rcomp = I0(:);
   RF_Tcomp = I1(:);
   sampleo=sps;
   save_tick=1;
if (save_tick==1)     
 files = ['RFs_' arsta{iar} '_' num2str(uorid(evt)) '.mat'];
 save  ZRT_data.mat  RF_Zcomp RF_Rcomp RF_Tcomp  istn iorid slowP baz  delta ...
        ymdtime arrtime dep MagB olat olon slat slon sampleo Ehang Nhang  % save data for receiver functions
if exist([savepath '/DatabySta_New1/' char(istn)],'dir')==0 
 unix(['mkdir ' savepath '/DatabySta_New1/' char(istn)]);
end
 unix(['mv ZRT_data.mat  ' savepath '/DatabySta_New1/' char(istn) '/' files]);  
 clear RF_Zcomp RF_Rcomp RF_Tcomp I0 I1 I2 istn iorid slowP baz delta ymdtime  ...
     arrtime dep MagB olat olon slat slon sampleo Ehang Nhang 
 trfree(tr); trdestroy(tr);
 end
 end % loop for each station 
%  if (save_tick==1)  
%  files = ['RFs_' num2str(uorid(evt)) '.mat'];
%  save  ZRT_data.mat   RF_Zcomp RF_Rcomp RF_Tcomp  istn iorid slowP baz  delta ...
%         ymdtime arrtime dep MagB olat olon slat slon sampleo Ehang Nhang  % save data for receiver functions
%  
%  unix(['mv ZRT_data.mat  ' savepath '/M5.5_New/' files]);  
%  clear RF_Zcomp RF_Rcomp RF_Tcomp I0 I1 I2 
%  trfree(tr);
%  end
fid=fopen([savepath '/DatabySta_New1/log.txt'],'a+');
fprintf(fid,'%d:%d \n',evt,uorid(evt));
fclose(fid);
end % loop for each event =================================================
% =========================================================================
dbclose(db);
