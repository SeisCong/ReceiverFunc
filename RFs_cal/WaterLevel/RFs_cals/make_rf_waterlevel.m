
% this program calculates ZRT receiver functions with waterlevel
% deconvolution

close all; clear all;
Workfolder = '/Users/congli/Project/ReceiverFunc_NewEngland/';
% datafolder = [Workfolder 'Data/RFs/'];
filename = dir([Workfolder 'Data/data_by_stations/*.mat']);

samprate = 20;
dt=1/samprate;
   
preA=-30;
postA=60;

fband=[0.05 2.0];
N=2;
Feff=1/dt/2;
fone=fband(1)/Feff; ftwo=fband(2)/Feff;
[b,a]=butter(N,[fone ftwo]);


! ls -d /Users/congli/Project/ReceiverFunc_NewEngland/Data/data_by_stations/* > StaList.txt

stalist = textread('StaList.txt','%s');

for ss = 3:length(stalist)
    
stadir = char(stalist{ss});   

sta = stadir(23:end)

%unix(['ls -d ' stadir '/M* > MagList.txt' ]);

%maglist = textread('MagList.txt','%s');

ff = 0;

unix([Workfolder 'Data/RFs/' sta]);

Z=[]; R=[]; T=[]; staname={}; baz=[]; orid=[]; slowP=[];
elat=[]; elon=[]; rlat=[]; rlon=[];
delta=[]; dep=[]; magB=[]; 

%for mm = 2:length(maglist)
    
filename = dir( [char(stalist{ss})  '/ZRT_*.mat'] );%lc

for i = 1:length(filename),
    
    ff = ff+1;

    clear files
    files = load( [char(stalist{ss}) '/' filename(i).name] );%lc
    
    R(:,ff) = files.RR;
    T(:,ff) = files.TT;
    Z(:,ff) = files.ZZ;
   
    staname{ff} = files.staname;

    baz(ff) = files.backazi;
    orid(ff) = files.orid;
    slowP(ff) = files.slowness;
    elat(ff) = files.elat;
    elon(ff) = files.elon;
    rlat(ff) = files.rlat;
    rlon(ff) = files.rlon;
    delta(ff) = files.dist;
    dep(ff) = files.depth;
    magB(ff) = files.magB;
    
end

%end  % of magnitude list



nevents = length(staname);     % number of stations for this event
npts = size(Z,1);
time = preA+(0:npts-1)*dt;

% define a tukey window
r=5; % secs
r=r/(dt*.5*npts);
Tf=tukeywin(npts,r);
    
Z0=zeros(size(Z)); R0=zeros(size(R)); T0=zeros(size(T));
  
% filter data
disp('filter the data ...')
for ff=1:nevents    
    Z0(:,ff)=filtfilt(b,a,Z(:,ff).*Tf);
    R0(:,ff)=filtfilt(b,a,R(:,ff).*Tf);
    T0(:,ff)=filtfilt(b,a,T(:,ff).*Tf);
end
offset=0;
% for ff=[33,18]    
% figure(1), hold on
% plot(Z(:,ff)./max(Z(:,ff))+offset,'k');
% % plot(Z0(:,ff)./max(Z0(:,ff))+offset,'r');
% offset=offset+5;
% end

%%% do water-level deconvolution
preA=30; 
Ns = (preA-25)*samprate;
Ne = (preA-5)*samprate;
time_shift = preA;
        
pa=time_shift*samprate;
       
% Slepian parameters, NW=P=main-lobe half-width (in freq bins)
NW  = 2.5;
nTap= 4;      % number of tapers used

disp('calculating receiver functions with waterlevel deconvolution ...')
for ff=1:nevents
    S=zeros(npts,1);

    S=Z(:,ff);
    RF_Z = dcv_revised( S, Z(:,ff), samprate, Ns, Ne, nTap,NW, time_shift );
    RF_R = dcv_revised( S, R(:,ff), samprate, Ns, Ne, nTap,NW, time_shift );
    RF_T = dcv_revised( S, T(:,ff), samprate, Ns, Ne, nTap,NW, time_shift );
    
    af = 2.5;
    RF_Zt = gaussf( RF_Z, samprate, af);
    RF_Rt = gaussf( RF_R, samprate, af);
    RF_Tt = gaussf( RF_T, samprate, af);
    
    Zraw=Z(:,ff); Rraw=R(:,ff); Traw=T(:,ff);
    Zfilt=Z0(:,ff); Rfilt=R0(:,ff); Tfilt=T0(:,ff);
    
    station=char(staname{ff});
    backazi=baz(ff);
    iorid=orid(ff);
    slowness=slowP(ff);
    eqlat=elat(ff);
    eqlon=elon(ff);
    stalat=rlat(ff);
    stalon=rlon(ff);
    dist=delta(ff);
    depth=dep(ff);
    mag=magB(ff);
    
    fname = ['RFs_waterlevel_' station '_' num2str(iorid) '.mat']
    fileex=exist(['data/RFs/' station '/'  fname],'file');
    if(fileex==2)
    eval(['save tmp_rf.mat Zraw Rraw Traw Zfilt Rfilt Tfilt RF_Z RF_R RF_T RF_Zt RF_Rt RF_Tt samprate station backazi iorid slowness eqlat eqlon stalat stalon dist depth mag'])
    unix(['mv tmp_rf.mat data/RFs/' station '/'  fname]);
    end
end
        
end % of stations        




