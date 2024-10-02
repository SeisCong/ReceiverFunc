% Interactive cross-correlation and plot for each event.

clear all; close all;

samprate=20;                % desired sample rate for data

workdir = './data_xcorr/';
mag = 'M5.4_5.5'
datadir = [workdir mag];
    
filename = dir([datadir '/*.mat']);
    
for fi = 1:length(filename), fi, filename(fi).name

    files = load( [datadir '/' filename(fi).name] ); 
    
    eorid = filename(fi).name(5:end-4);
    
    Z = files.Zx_raw;
    R = files.Rx_raw;
    T = files.Tx_raw;
    iorid = files.iorid_new;
    istn = files.istn_new;
    baz = files.baz_new;
    dep = files.dep_new;
    slowP = files.slowP_new;
    MagB = files.MagB_new;
    delta = files.delta_new;
    ymdtime = files.ymdtime_new;
    slat = files.slat_new;
    slon = files.slon_new;
    olat = files.olat_new;
    olon= files.olon_new;
    
    nstn = length(istn);
    
    for nn = 1:nstn
     files = ['ZRT_' char(istn{nn}) '_' num2str(eorid) '.mat'];
     
     ZZ=Z(:,nn);
     RR=R(:,nn);
     TT=T(:,nn);
     staname=char(istn{nn}); 
     orid=iorid(nn);
     slowness=slowP(nn);
     backazi=baz(nn);
     dist=delta(nn); 
     etime=ymdtime{nn};
     depth=dep(nn);
     magB=MagB(nn);
     elat=olat(nn);
     elon=olon(nn);
     rlat=slat(nn);
     rlon=slon(nn); 
     
     save  tmp_data.mat   ZZ RR TT  staname orid slowness backazi  dist ...
        etime depth magB elat elon rlat rlon samprate  % save data for receiver functions

    targetpath = ['./data_by_stations/' char(istn{nn}) '/' mag]; 
    if (exist(targetpath, 'dir') == 0), mkdir(targetpath); end
    
    unix(['mv tmp_data.mat ' targetpath '/'  files]);  
   end
        
 
end % loop for each event 
