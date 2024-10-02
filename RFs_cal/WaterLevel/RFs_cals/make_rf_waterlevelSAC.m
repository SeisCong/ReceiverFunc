%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE IS USED FOR EXTRACTING DATA FROM SAC FILES
% AND THEN MAKE DECONVOLUTION 
%
%
%  CONG LI
% 06/16/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear;close all;clc;
 %run('/opt/antelope/5.5/setup.m');
 %! ls -d /Users/congli/RFs_Research/RFs_NewEngland/NewData/* > StaListNewData.txt;
 stalist = textread(['StaListNewData.txt'],'%s');
 nstn = length(stalist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filter Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
samprate = 20;
dt=1/samprate;
fband=[0.1 1.5];
N=2;
Feff=1/dt/2;
fone=fband(1)/Feff; ftwo=fband(2)/Feff;
[b,a]=butter(N,[fone ftwo]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Deconvolution Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PreA=30;
PostA=100;
Ns = (PreA-25)*samprate;
Ne = (PreA-5)*samprate;
time_shift = PreA;        
pa=time_shift*samprate;
NW  = 2.5;    % Slepian parameters, NW=P=main-lobe half-width (in freq bins)
nTap= 4;      % number of tapers used
 for i=12:length(stalist) 
 disp('new station');
 fnamez=dir([stalist{i} '/*.z']);
 for j=1:length(fnamez(:,1))
 filenamez=[stalist{i} '/' fnamez(j).name(1:22) '.z'];
 filenamer=[stalist{i} '/' fnamez(j).name(1:22) '.r'];
 filenamet=[stalist{i} '/' fnamez(j).name(1:22) '.t'];
 outputz = rsac(filenamez);
 outputr = rsac(filenamer);
 outputt = rsac(filenamet);
 if isempty(outputz)||isempty(outputr)||isempty(outputt)
     disp([fnamez(j).name(1:22) 'does not have 3 components'])
     continue;
 end     
 tt=outputz(:,1);
 Zraw=outputz(:,2);
 Rraw=outputr(:,2);
 Traw=outputt(:,2);
 header=outputz(:,3);
 backazi=header(53);
 spst=header(1);
 sps=round(1/spst);
 dist=header(54);
 distkm=header(51);
 stalat=header(32);
 stalon=header(33);
 eqlat=header(36);
 eqlon=header(37);
 depth=header(39)/1000;
 mag=header(40);
 npts=header(80);
 P_arr=header(11);
 station1=char(fnamez(j).name(19:22));
 station=['X8.' station1];
 date=fnamez(j).name(1:11);
 if (P_arr<=0)
     continue;
 end
 tauP=['taup_time -mod iasp91 -h ' num2str(depth) ' -ph P -deg ' num2str(dist)];
 unix(tauP);
AntelopePath   % set Antelope path
run('/opt/antelope/5.5/setup.m');
slowJ=[];phasenames=[];
[slowJ, phasenames] = arr_slowness(dist,depth);
numPhase=length(slowJ);
for phaseE=1:numPhase       
      pp=cell2mat(phasenames(phaseE));
      if length(pp)==1 && pp(1)=='P'
        phaseNum=phaseE;  
      end
end  
slowness=slowJ(phaseNum);  % assumes phase is P or PKIKP
if (dist>30&&dist<95)
%%%%%%%%%%%%%%%%Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Filter the data ...')
% define a tukey window
r=5; % secs
r=r/(dt*.5*npts);
Tf=tukeywin(npts,r);
Zfil=zeros(size(Zraw)); Rfil=zeros(size(Rraw)); Tfil=zeros(size(Traw));
Zfil=filtfilt(b,a,Zraw.*Tf);
Rfil=filtfilt(b,a,Rraw.*Tf);
Tfil=filtfilt(b,a,Traw.*Tf);
%%%%%%%%%%%%%%%%Cut%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zcut=Zfil(round((P_arr-tt(1)-PreA)*sps) : round((P_arr-tt(1)+PostA)*sps));
Rcut=Rfil(round((P_arr-tt(1)-PreA)*sps) : round((P_arr-tt(1)+PostA)*sps));
Tcut=Traw(round((P_arr-tt(1)-PreA)*sps) : round((P_arr-tt(1)+PostA)*sps));
%%%%%%%%%%%%%%%%Resample%%%%%%%%%%%%%%%%%%%%%%%%%
Zres = resample(Zcut,samprate, sps);
Rres = resample(Rcut,samprate, sps);
Tres = resample(Tcut,samprate, sps);
npts1=length(Zres);
S=zeros(npts1,1);
S=Zres;
RF_Z = dcv_revised( S, Zres, samprate, Ns, Ne, nTap,NW, time_shift );
RF_R = dcv_revised( S, Rres, samprate, Ns, Ne, nTap,NW, time_shift );
af = 2.5;
RF_Zt = gaussf( RF_Z, samprate, af);
RF_Rt = gaussf( RF_R, samprate, af);
unix(['mkdir data/RFs/NewData/' station]);
fname = ['RFs_waterlevel_' station '_' date '.mat'];     
eval(['save tmp_rf.mat Zraw Rraw Traw Zres Rres Tres RF_Z RF_R RF_Zt RF_Rt samprate tt station backazi date eqlat eqlon stalat stalon dist slowness depth mag'])
unix(['mv tmp_rf.mat data/RFs/NewData/' station '/' fname]); 
end
end
aaaaa=1;
end