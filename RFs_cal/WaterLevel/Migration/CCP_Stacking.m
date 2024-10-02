%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CCP STACKING CODE
% The common conversion point (CCP) stacking of the RFs in the depth domain
% YOU CAN CHANGE THE SMOOTHNESS OF YOUR IMAGE DURING THE STACKING BY
% ADJUSTING THE BIN SPACING AND THE RADIUS OF THE BINS THAT YOU USE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modified history
% 3/12/2016 Haiying Gao, 
% 7/21/2017 Cong Li, Make the Iasp91 as a defualt model, imported model as a updated model
% 
% Contact:Cong Li,conli@geo.umass.edu
% Umass Amhers, MA 01003, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all
addpath('/Users/congli/Software/Mat_Lib/');
addpath('/Users/congli/Software/Mat_Lib/CCP');
workfolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/';
RFsfolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/DataBySta_New_autoselect/';
Resultfolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/CCP/';
Pickingfolder = [workfolder 'Picking/'];
statname = textread([workfolder 'StaList_NewData.txt'],'%s'); % Read the staion list
%[statname,Stalat,Stalon,t1,t2,t3]=textread([Pickingfolder 'SingleTracePick_final.txt'],'%s %f %f %f %f %f');
hh = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('==========stack_moho===================== \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxpt = 5000;% max sampling number, may be changed based on the data
% samprate=20;
tsca = 0; 
DtoR = pi/180;
 icomp  = 'R'    ;% Component flag
 imove = 0       ;% moveout flag (0=Pds; 1=2p1s; 2=2s1p)
 numitera = 3    ;% Bootstrap iteration,use 200 for final
 Ampcut = 0.000   ;% threshold value for recording the max amplitude
 PWOrder = 0.0    ;% weighting filter power (0=linear; <1=gentler; >1=stronger)
 wl   = 0.05      ;% low pass corner
 wh   = 2.0      ;% hi pass corner
 wbaz1 = 0       ;% lower baz for windowing 90/180/180:260/270:360
 wbaz2 = 360     ;% upper baz for windowing
 dzi = 1         ;% depth increment, use 1 for final run
 dzmax = 150.0   ;% max depth
 rotang = 0*DtoR ;%stack into a set of bins that are not oriented n-s and e-w
 refAng = 20     ;%reference P-angle of incidence for amplitude scaling
 reModel= 1     ; %reference Model, if yes(1), need model input, otherwise it will use iasp91 model as the defualt model
 guass_af=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  LOAD DATA                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh = 1;
prf=[];seis=[]; tseis=[]; istn={}; slat=[]; slon=[];
baz=[];slow=[]; olat=[];  olon=[]; depth=[]; Rmt=[]; Tmt=[];
time=[]; ttime=[];iorid=[];
for ss=1:length(statname)
datafolder = [char(statname{ss}) '/'];
filelist= dir([RFsfolder datafolder '*.mat']);
for i = 1:length(filelist)  
   clear files 
   files = load( [RFsfolder datafolder filelist(i).name] );
   if (length(files.RF_Zcomp)< maxpt)
       max_pt=length(files.RF_Zcomp);   
       pseis(1:max_pt,hh) = files.RF_Zcomp(1:max_pt);
       sseis(1:max_pt,hh) = files.RF_Rcomp(1:max_pt);
       pseis(max_pt+1:maxpt,hh)=0;
       sseis(max_pt+1:maxpt,hh)=0;
   else
       pseis(:,hh) = files.RF_Zcomp(1:maxpt);
       sseis(:,hh) = files.RF_Rcomp(1:maxpt);
   end
   iistn{hh} = files.istn;
   mpt = find(pseis(:,hh)==max(pseis(:,hh)));
   npt= length(pseis(:,hh));
   samprate= files.sampleo;
   time(:,hh) = (0:npt-1)/samprate - mpt(1)/samprate;  
   bbaz(hh)   = files.baz;
   sslow(hh)  = files.slowP;
   sslat(hh)  = files.slat;
   sslon(hh)  = files.slon;
   oolat(hh)  = files.olat;
   oolon(hh)  = files.olon;
   iiorid(hh)  = files.iorid;
   hh=hh+1;  
end
end
clear files
[npts,ntra] = size(sseis);
dt=1/samprate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    GROUP DATA IN TERMS OF BACK-AZIMUTH, DECREASE THE AZIMUTH EFFECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = find(bbaz>=wbaz1 & bbaz<=wbaz2);
iorid=zeros(1,length(index));
for kk=1:length(index)
 t0=find(time(:,kk)==0);
tlength = 2000;
if t0 >= (size(sseis,1)-tlength), continue, end
npre =t0 ; npost =t0+tlength; % start from P arrivals and 50secs after P
if sseis(t0,index(kk))<0  %% Whether the negative polarity needs to be corrected
%     pseis(:,index(kk))=-pseis(:,index(kk));
%     sseis(:,index(kk))=-sseis(:,index(kk));
end
prf(:,kk) = pseis(npre:npost,index(kk));
seis(:,kk) = sseis(npre:npost,index(kk));
istn{kk} = iistn(index(kk));
slat(kk) = sslat(index(kk));
slon(kk) = sslon(index(kk));
baz(kk) = bbaz(index(kk));
slow(kk) = sslow(index(kk));
olat(kk) = oolat(index(kk));
olon(kk) = oolon(index(kk));
if (exist('iiorid','var'))
% iorid(kk)= iiorid(index(kk));
end
prf(:,kk) = gaussf(prf(:,kk), samprate, guass_af);
seis(:,kk) = gaussf(seis(:,kk), samprate, guass_af);
end
clear pseis sseis iistn sslat sslon bbaz sslow oolat oolon ttime iiorid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      FILTER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b,a,Tf]=GenerateFilter(seis,dt,wl,wh);
for i = 1:index
 prf(:,i)=filtfilt(b,a,prf(:,i)).*Tf;
 seis(:,i)=filtfilt(b,a,seis(:,i)).*Tf;
 seis(:,i)  = seis(:,i)/max(prf(:,i));%max(seis(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Grid Parameter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dZ = 0:dzi:dzmax;% depth vector
sdcut = 0.30;% max std. dev based on trace stdev
supercut = 1;% min # traces to keep bin 
dmx = 20;%x bin center spacing, control spacing between bin center points in the x and y directions
dmy = 20;%y bin center spacing
binsz1 = mean([dmy dmx])*1.5;% bin radius as a function of bin center spacing
binsz2 = mean([dmy dmx])*3.0;% bin radius as a function of bin center spacing
binnum_x = 90; binnum_y =90; % number of bins
xmax = (binnum_x/2)*dmx; xmin = -(binnum_x/2)*dmx;
ymax = (binnum_y/2)*dmy; ymin = -(binnum_y/2)*dmy;
Rv = 1.76; % Vp/Vs value
ax   = -76;% longitude of array center
ay   = 40.0; % latitude of array center
xsta = (slon-ax)*111.13.*cos(ay*pi/180); ysta = (slat-ay)*111.13; % location(distance from the center (0,0))
if ntra >= 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                CCP STACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Remap RF time to vertical incidence time (ie., RF migration) \n')
% Time to Depth transformation
ndz = size(dZ,2); zthk = ones(1,ndz)*dzi;
[nseis,svel,mean_D] = DepthTrans(dZ,dzi,dt,reModel,ntra,slon,slat,slow,seis,imove);
clear seis;
% Generate grids
fprintf('CALC x,y,z bin crds \n')
ymod = [ymin:dmy:ymax]; % Grids in y direction
xmod = [xmin:dmx:xmax]; % Grids in x direction
nx = size(xmod,2); ny = size(ymod,2);
iy = 0; 
for yy = ymod% Mesh 
  iy = iy  + 1; ix = 0;
  for xx = xmod
    ix = ix  + 1;
    bin(iy,ix,1) = xx;
    bin(iy,ix,2) = yy;
  end
end
% Calculate the stacked RFs index and weight in each grid
mean_D=mean_D*1.5;
[wt,nnxy,inxy,scale] = CalWeight(xmod,ymod,bin,xsta,ysta,ntra,dZ,dzi,...
refAng,DtoR,svel,slow,baz,mean_D,supercut);
nseis=nseis.*scale;% Amp scale correction to the designed angle
% CMP Stacking
[cmp,sdev,AMPmax,DEPmax] = CMPStacking( dzmax,dzi,dZ,numitera,nx,ny,...
nnxy,inxy,wt,nseis,supercut,Ampcut,PWOrder);
% make and save output results file name 
ofile = ['CCP_US2016_' num2str(dmx) '_'  num2str(dzmax) '_' icomp...
    '.moho_' num2str(imove) '_' num2str(360+ax) ...
    '_' num2str(ay) '_' num2str(wl) '_' num2str(wh) '_' num2str(PWOrder) ...
    '_0-360' 'testRFs.mat']
fprintf('OUTPUT file= %s \n',ofile)

eval(['save ' Resultfolder  ofile ' cmp sdev ymod xmod dZ ndz icomp wl wh sdcut dmy dzi ' ...
  ' wbaz1 wbaz2 ntra nx ny imove nnxy binsz1 ' ...
  'ax ay xsta ysta istn bin xmax xmin ymax ymin dmx dmy inxy nseis baz'])
end
clear all;