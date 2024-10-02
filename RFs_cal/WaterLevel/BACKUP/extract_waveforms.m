% Interactive cross-correlation and plot for each event.

clear all; close all;
AntelopePath   % set Antelope path

samprate=20;                % desired sample rate for data

%===== extract data from Antelope database =====

 network = {'RFs_1995_now'};

 dir = ['/Users/congli/Research/RFs_NewEngland/' char(network) '/' char(network) ]
 db = dbopen( dir, 'r' );
 dbar = dblookup_table( db, 'arrival' );
 dbar = dbsubset( dbar, 'iphase != "del"');
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
 [artime,arsta,bazimuth,azimuth,arid,orid,edep,lat_ori,lon_ori,del,Mb] = ...
    dbgetv(dbar,'time','sta','seaz','esaz','arid','orid','depth','lat','lon','delta','mb');
%  strtime(artime)
uorid = unique(orid);
ind = 1;
display(['number of events ' num2str(length(uorid))]);

Phase=zeros(length(uorid),1);
for evt = 1:length(uorid) % Find Mean slowness of P or PKIKP wave based on iasp91 model
%for evt = 1
%  index = find(orid == uorid(evt));
  index11 = find(orid == uorid(evt));  % changed from the origonal program get rid of repeated stations for each event
  [trash, index22] = unique(arsta(index11));
  index = index11(index22); 
  
  bazA(evt)=mean(bazimuth(index));
  delA(evt)=mean(del(index));
  depA(evt)=mean(edep(index));
  [slowJ, phasenames] = arr_slowness(delA(evt),depA(evt));% <----
  numPhase=length(slowJ);
  for phaseE=1:numPhase
    if Phase(evt)==0 % find which occurs first, P or PKPdf
        
      pp=cell2mat(phasenames(phaseE));
      if length(pp)==5 & pp(1:5)=='PKPdf';
        phaseNum=phaseE;
        Phase(evt)=2;
      elseif length(pp)==1 && pp(1)=='P'
        phaseNum=phaseE;
        Phase(evt)=1;
      end
      
    end
  end
  
  slowA(evt)=slowJ(phaseNum);  % assumes phase is P or PKIKP
end
slowA=slowA*111.195;

preA=-30;
postA=100;
nt=(postA-preA)*samprate; % number of time points  
   
% =========================================================================
for evt = 1:length(uorid)
%for evt = 1
   fprintf('Origin Id %d: %d\n',evt, uorid(evt));
      
   R=[]; T=[]; Z =[]; slowP=[]; sta_coor=[]; ymdtime={}; istn={};
   iorid=[]; olat=[]; olon=[]; slat=[]; slon=[]; baz=[]; atime=[]; dist=[];
   arrivaltime=[]; dep=[];Rx=[];Tx=[];Zx=[];MagB=[];
   
   index11 = find(orid == uorid(evt));  % changed from the origonal program get rid of repeated stations for each event
   [trash, index22] = unique(arsta(index11));
   index = index11(index22);
   nStn = length(index);     % number of stations for this event

   count = 0;
   for ind =1:nStn   %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   iar = index(ind);
   [slowness, phasenames] = arr_slowness( del(iar), edep(iar));
   slowP(ind) = slowness(1,1);
   
   dbsub = dbsubset( dbwf, ['sta == "' arsta{iar} '"'] );
   tr = trload_css(dbsub, artime(iar)+preA, artime(iar)+postA );% <---
   trsplice( tr, 0.5 );
   trapply_calib(tr);
   nrecs = dbnrecs(tr);

    if ( dbnrecs(tr) ~= 3 )
        disp('Not 3 components ...'); 
        continue
    else
        tr.record = 0;
        H1 = trextract_data(tr);
        tr.record = 1;
        H2 = trextract_data(tr);
        tr.record = 2;
        HZ = trextract_data(tr);
        sps = dbgetv( tr, 'samprate' );
        if sps~=round(sps)
            sps=round(sps);
        end

        if length(H1)~=length(H2), continue, end
        if ~isempty(find(isnan(HZ))), continue, end
           
        HZ = resample(HZ,samprate, sps);
        H1 = resample(H1,samprate, sps);
        H2 = resample(H2,samprate, sps);
        HZ = detrend(HZ);
        H1 = detrend(H1);
        H2 = detrend(H2);
    end
  
   nrecs = dbnrecs(tr);
   dbsub = dbsubset( dbst, ['sta == "' arsta{iar} '"'] );
   [slt,sln,elv]  = dbgetv( dbsub, 'lat','lon','elev');
   
   dbsub = dbsubset( dbnw, ['sta == "' arsta{iar} '"'] );
   nwsta  = dbgetv( dbsub, 'snet');
   
   dbsub = dbsubset( dbchan, ['sta == "' arsta{iar} '" && ' sprintf(' orid ==%f ', uorid(evt)) '']);
   %dbsub = dbsubset( dbchan, 'sta == "D61A" && orid == 21498');
   [hang,staname] = dbgetv( dbsub, 'hang','sta');
   if (hang(1)~=90)
       ;
   end
   if isempty(char(arsta(iar))), continue, end
   
   if length(HZ)~=2600, continue, end
   if length(H1)~=2600, continue, end
   
   count = count+1;
   iorid(count) = orid(iar);
   istn(count,:)= {[char(nwsta) '.' char(arsta(iar))]};
   olat(count)  = lat_ori(iar);
   olon(count)  = lon_ori(iar);
   slat(count)  = mean(slt);
   slon(count)  = mean(sln);
   selv(count)  = mean(elv);
 if (size(hang,1) == 3)
     Ehang(count)=hang(1);
     Nhang(count)=hang(2);
     Zhang(count)=hang(3);
 elseif (size(hang,1)==2)
     Ehang(count)=hang(1);
     Nhang(count)=hang(2);
     Zhang(count)=9999;    % 9999 means unreliable value
 elseif (size(hang,1)==1)
     Ehang(count)=hang(1);
     Nhang(count)=9999;
     Zhang(count)=9999;
 else
     Ehang(count)=9999;
     Nhang(count)=9999;
     Zhang(count)=9999;
 end
%    hangE(count)=hang(1);
%    selv_corr(count) = round( (selv(count)/6.5)*samprate );
   baz(count)   = bazimuth(iar);
   delta(count)  = del(iar);
   dep(count) = edep(iar);
   MagB(count) = Mb(iar);
   atime(count) = artime(iar);
   ymdtime{count} = strydtime(atime(count)); % convert arrival time to hhmmss
   aa = strydtime(atime(count)); % convert arrival time to hhmmss
   arrt = aa(18:end);
   arrtime(count) = str2num(arrt(1:2))*3600 + ...
                      str2num(arrt(4:5))*60 + ...
                      str2num(arrt(7:8)) + ...
                      str2num(arrt(10:end))*0.001;
   Zcomp(:,count) = HZ(:);
   H1comp(:,count) = H1(:);
   H2comp(:,count) = H2(:);
   
   end % loop for each station 
     
 files = ['ZRT_' num2str(uorid(evt)) '.mat'];
 save  ZRT_data.mat   Zcomp H1comp H2comp  istn iorid slowP baz  delta ...
        ymdtime arrtime dep MagB olat olon slat slon samprate Ehang Nhang  % save data for receiver functions

 unix(['mv ZRT_data.mat data/M5.0/'  files]);  

end % loop for each event =================================================
% =========================================================================

dbclose(db);
