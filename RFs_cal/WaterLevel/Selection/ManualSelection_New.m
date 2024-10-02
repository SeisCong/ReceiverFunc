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
addpath('/Users/congli/Software/Mat_Lib');
workfolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/'
%datafolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/DataBySta_New_autoselect/'
datafolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/DataBySta_New/'
Resultfolder ='/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/Picking/'
stalist = textread(['../StaList_NewData.txt'],'%s'); % Read the staion list
nstn = length(stalist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%samprate=40; % sampling rate
guass_af=10;
maxpt = 5000; % npts
ReferModel_tick=0;% move-out reference model tick
Pick_tick=1;% Picking tick
QuadrantStack_tick=1; % QuadrantStack tick
MeanStack_tick=1; % MeanStack tick
Moveout_tick=1;% moveout tick
PmsWindow=[2.5,6.5]; % Window for picking Pms
%fband = [0.05,0.75]; % Window for filter, try different frequencies
%fband = [0.1,1.2];
fband = [0.1,1.2];
%fband = [0.2,1.2];
%fband = [0.05,2.0];

uncertainty=1.0; % maximum limit for Pms picking arrival
Show_num=50;
amp=3;
del_num=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss =154:nstn, ss
station = char(stalist(ss));
redraw=0;
while(redraw==0) % redraw
filename = dir( [datafolder station '/RFs_*.mat'] );
RFs_Z=[]; RFs_R=[]; staname={}; baz=[]; iorid=[];stalat=[];stalon=[];slowness=[];RFs_Zt=[];
RFs_Rt=[];RFs_Tt=[];eqlat=[];eqlon=[];depth=[];mag=[];dist=[];
k=1; 
%del=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the RFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(filename)
   files = load( [datafolder station '/' filename(i).name] );
   if length(files.RF_Zcomp) < maxpt
       continue;
   end
   RFs_Zt(i,:) = files.RF_Zcomp(1:maxpt);
   RFs_Rt(i,:) = files.RF_Tcomp(1:maxpt);
   RFs_Tt(i,:) = files.RF_Rcomp(1:maxpt); 
   staname{i} = files.istn; 
   baz(i) = files.baz;
   iorid(i) = files.iorid;
   stalat(i)=files.slat;
   stalon(i)=files.slon;
   slowness(i)=files.slowP;
   samprate = files.sampleo;
   eqlat(i)=files.olat;
   eqlon(i)=files.olon;
   depth(i)=files.dep;
   mag(i)=files.MagB;
   dist(i)=deg2km(distance(eqlat(i), eqlon(i), stalat(i), stalon(i)));
end
dt=1/samprate;
RFs_Z = gaussf( RFs_Zt', samprate, guass_af);
RFs_R = gaussf( RFs_Rt', samprate, guass_af);
RFs_Z=RFs_Z';
RFs_R=RFs_R';
numfile=size(RFs_Z,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Time Move-out Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Moveout_tick~=1
RFs_Rm=RFs_R;
RFs_Zm=RFs_Z;
elseif Moveout_tick==1
 RFs_Rm=[];  
 RFs_Zm=[];
for i=1:numfile
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
set(gca,'XLim',[-2 15])
set(gca,'YLim',[cc-1 cc+Show_num])
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
title([ '(c) Radial RFs for Station ' char(staname{2})],'FontSize',20 )
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
if max(tmpZm)==0
    continue;
end
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
   tt= find( time >= 30 );  
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
   for p_max=1:50 
   next_del=input('\n next delet (Yes,0)?')';
   if next_del==0
       [del_x,del_y]=ginput(1);
       del_index=round(del_y)
       for p_del=1:length(time)
           del_plot(p_del)=del_index;
       end
       plot(time,del_plot,'b','LineWidth',2);
       del_num=del_num+1;
        del(del_num)=del_index;
    else
        break;
    end
   end
   close all;
%    pause;
end % end of loop for each event

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Manually select RFs and Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
next=input('\n Next Station (Yes,1/No,0)?')' % uncomment if you want to select;
if next==0
% del=input('\n\nDelete Trace Number:');
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
%     clear all;
    continue;    
end
del=[];del_num=0;
close all
end
end