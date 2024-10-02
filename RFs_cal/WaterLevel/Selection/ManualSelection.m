%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort RFs by the back-azimuth, manually select the RFs 
%
%% Cong Li
% History
% 05/16/2016
% 06/03/2017 
%       Time Move-out correction, based on a iasp91 model 
% 10/11/2017
%       Merge Sac data into the antelope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all
addpath('/Users/congli/Software/Mat_Lib/Moveout');
workfolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/'
datafolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/DataBySta_New/'
Resultfolder ='/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/Picking/'
stalist = textread(['../StaList.txt'],'%s'); % Read the staion list
nstn = length(stalist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samprate=40; % sampling rate
guass_af=10;
maxpt = 2000; % npts
ReferModel_tick=0;% move-out reference model tick
Pick_tick=1;% Picking tick
QuadrantStack_tick=1; % QuadrantStack tick
MeanStack_tick=1; % MeanStack tick
Moveout_tick=1;% moveout tick
PmsWindow=[2.5,6.5]; % Window for picking Pms
%fband = [0.05,0.75]; % Window for filter, try different frequencies
fband = [0.1,1.2];
%fband = [0.2,2.0];
uncertainty=1.0; % maximum limit for Pms picking arrival
Show_num=50;
amp=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss =1:nstn, ss
station = char(stalist(ss));
redraw=0;
while(redraw==0) % redraw
filename = dir( [datafolder station '/RFs_*.mat'] );
RFs_Z=[]; RFs_R=[]; staname={}; baz=[]; iorid=[];stalat=[];stalon=[];slowness=[];
k=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the RFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(filename)
   files = load( [datafolder station '/' filename(i).name] );
   if length(files.RFsZ_raw{1}) < maxpt
       continue;
   end
   RFs_Zt(i,:) = files.RFsZ_raw{1}(1:maxpt);
   RFs_Rt(i,:) = files.RFsR_raw{1}(1:maxpt);
   RFs_Tt(i,:) = files.RFsT_raw{1}(1:maxpt); 
   staname{i} = files.staname; 
   baz(i) = files.backazi;
   iorid(i) = files.orid;
   stalat(i)=files.rlat;
   stalon(i)=files.rlon;
   slowness(i)=files.slowness;
   samprate = files.samprate;
   eqlat(i)=files.elat;
   eqlon(i)=files.elon;
   depth(i)=files.depth;
   mag(i)=files.magB;
end
dt=1/samprate;
RFs_Z = gaussf( RFs_Zt, samprate, guass_af);
RFs_R = gaussf( RFs_Rt, samprate, guass_af);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Time Move-out Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Moveout_tick~=1
RFs_Rm=RFs_R;
RFs_Zm=RFs_Z;
elseif Moveout_tick==1
 RFs_Rm=[];  
 RFs_Zm=[];
for i=1:length(filename)
 [RFs_Rm(i,:),RFs_Zm(i,:)]=TimeMoveOut_Iasp91(RFs_R(i,:),RFs_Z(i,:),slowness(i),dt,maxpt);
end   
end
npts=size(RFs_Z,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wn=fband.*(dt*2); 
[b,a]=butter(2,Wn);  % create filter coefficients (second order)
r=0.5; % secs
r=r/(dt*.5*npts);
Tf=tukeywin(npts,r);
stackedR=0;stackedRm=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Ploting & Selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nevents = length(baz);
trace=nevents;
[ixx,baz_ind]  = sort(baz);
for cc=1:Show_num:nevents
figure('Position',[200 400 800 1600]), hold on, box on, grid on
set(gca,'FontSize',15);
xlabel('Time after P-wave (s)', 'FontSize',15)
set(gca,'Xtick',[-10:1:20],'ygrid','on')
set(gca,'XLim',[-2 10])
set(gca,'YLim',[cc-1 cc+Show_num])
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
title([ '(c) Radial RFs for Station ' char(staname{1})],'FontSize',20 )
for j = cc:cc+Show_num-1    
   if (j>trace)
       break;
   end
   tmpZm=[]; tmpRm=[];
   tmpZm=RFs_Zm(baz_ind(j),:); % move-out RFs
   tmpRm=RFs_Rm(baz_ind(j),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter and normalization        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   tmpZm=filtfilt(b,a,tmpZm').*Tf;
   tmpRm=filtfilt(b,a,tmpRm').*Tf;  
   tmpRm=amp*tmpRm/max(abs(tmpRm));
   tmpZm=tmpZm/max(abs(tmpZm));   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get time range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mpt = find(tmpZm==max(tmpZm));
   expect_mpt=samprate*30;
   if mpt~=expect_mpt
       disp(['mpt != ' num2str(expect_mpt) ])
   end
   time = (0:npts-1)/samprate - mpt/samprate;
   tt= find( time >= 10 );  
    if (tmpRm(mpt)>0)
        ;
    else
        tmpRm=-tmpRm;
    end
    RFs_Rm_pick(j,:)=tmpRm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Drawing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   rr1m = tmpRm +j; rr1m(tt) = j;
   rr1m(1)=j; rr1m(end)=j;
   fii1m = find(rr1m <= j); frr1m = rr1m; frr1m(fii1m) = j;
   fii2m = find(rr1m >= j); frr2m = rr1m; frr2m(fii2m) = j;
   plot(time,rr1m,'k','LineWidth',1); 
   fill(time,frr1m','r','facealpha',0.7);
   text(-3.5, j+0.2, [num2str(round(baz(baz_ind(j))))],'FontSize',15)
   text(-2.5, j+0.2, [num2str(j)],'FontSize',15)
end
   drawnow
   pause;
end % end of loop for each event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Manually select RFs and Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
next=input('\n Next Station (Yes,1/No,0)?')' % uncomment if you want to select;
if next==0
del=input('\n\nDelete Trace Number:');
   if (length(del)~=0)
    temptpath = [datafolder station '/' 'tmp1']; 
    if (exist(temptpath, 'dir') == 0)
        mkdir(temptpath);
    end
    for i=1:length(del)
     unix(['mv ' datafolder station '/' filename(baz_ind(del(i))).name ' ' temptpath]);
    end
    end
else
    redraw=1;
    close all;
    continue;
end
close all
end
end