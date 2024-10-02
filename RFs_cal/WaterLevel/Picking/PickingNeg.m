%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Picking the arrival time and the range for Pms,stack RFs based on the Quandrant
%% Cong Li
% History
% 05/16/2016
% 06/03/2017 
%       Time Move-out correction, based on a iasp91 model 
% 10/11/2017
%       Merge Sac data into the antelope
%       Pick the each selected RFs as the uncertainty
%       Stacked the RFs based on the back-azimuth 
% 10/30/2017
%       Pick the stacked PmS arrival
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all
addpath('/Users/congli/Software/Mat_Lib/Moveout');
addpath('/Users/congli/Software/Mat_Lib/Stack');
workfolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/';
datafolder = '/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/DataBySta_New_autoselect/';
Resultfolder ='/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/Picking/';
% Savepath ='/Users/congli/Project/ReceiverFunc_Appalachians/Database/Data/'
stalist = textread([workfolder 'StaList_NewData.txt'],'%s'); % Read the staion list
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
PmsWindow=[2.5,4.0]; % Window for picking Pms
%PmsWindow=[5.0,9.0]; % Window for picking Pms
%fband = [0.05,0.75]; % Window for filter, try different frequencies
fband = [0.05,1.2];
%fband = [0.2,2.0];
uncertainty=1.0; % maximum limit for Pms picking arrival
azi=1;
amp=20;
fid = fopen([Resultfolder 'SingleTracePick_neg_Inner.txt'], 'a+'); % Record of Pms arrival and range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss =538%:nstn, ss
station = char(stalist(ss));
redraw=0;
while(redraw==0) % redraw
filename = dir( [datafolder station '/RFs_*.mat'] );
if length(filename)<=1
    break;
end
RFs_Z=[]; RFs_R=[]; staname={}; baz=[]; iorid=[];stalat=[];stalon=[];slowness=[];
k=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the RFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(filename)
   files = load( [datafolder station '/' filename(i).name] );
   if length(files.RF_Zcomp) < maxpt
       continue;
   end
   RFs_Zt(i,:) = files.RF_Zcomp(1:maxpt);
   RFs_Rt(i,:) = files.RF_Rcomp(1:maxpt);
   RFs_Tt(i,:) = files.RF_Tcomp(1:maxpt); 
   
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
   
end
dt=1/samprate;
RFs_Z = gaussf( RFs_Zt', samprate, guass_af);
RFs_R = gaussf( RFs_Rt', samprate, guass_af);
RFs_Z=RFs_Z';
RFs_R=RFs_R';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Synchronize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expect_mpt=samprate*30+1;
for i=1:size(RFs_R,1)
   mptZ = find(RFs_Z(i,:)==max(RFs_Z(i,:)));
   movePZ=mptZ-expect_mpt;
   if(movePZ~=0)
       if (movePZ<0)
           RFs_Z(i,:)=[zeros(1,abs(movePZ)),RFs_Z(i,1:end+movePZ)];
       elseif (movePZ>0)
           RFs_Z(i,:)=[RFs_Z(i,movePZ+1:end),zeros(1,abs(movePZ))];
       end
   end
   mptR = find(RFs_R(i,:)==max(RFs_R(i,:)));
   movePR=mptR-expect_mpt;
   if(movePR~=0)
       if (movePR<0)
           RFs_R(i,:)=[zeros(1,abs(movePR)),RFs_R(i,1:end+movePR)];
       elseif (movePR>0)
           RFs_R(i,:)=[RFs_R(i,movePR+1:end),zeros(1,abs(movePR))];
       end   
   end
end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stack and mean within a certain degree, reduce the uncertainty of picking time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RFs_Ram=[];RFs_Zam=[];Tpick=[];
nevents = length(baz);
%  [ixx,baz_ind]  = sort(baz);
%  [bazShow,bazShowind]=unique(azi*floor(baz(baz_ind)/azi));
if MeanStack_tick==1;
[ RFs_Ram,RFs_Zam,baz_ind,bazShow,bazShowind ] = Stack_BazMean(RFs_Rm,RFs_Zm,baz,maxpt,azi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wn=fband.*(dt*2); 
[b,a]=butter(2,Wn);  % create filter coefficients (second order)
r=0.5; % secs
r=r/(dt*.5*npts);
Tf=tukeywin(npts,r);
stackedR=0;stackedRm=0;
qua1_num=0; qua2_num=0; qua3_num=0; qua4_num=0;
stackedRm_Q1=0; stackedRm_Q2=0; stackedRm_Q3=0; stackedRm_Q4=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Ploting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlength=30;
figure('Position',[200 400 800 1600]), hold on, box on, grid on
set(gca,'FontSize',10);
xlabel('Time after P-wave (s)', 'FontSize',15)
set(gca,'Xtick',[-10:5:30],'ygrid','on')
set(gca,'XLim',[0 15])
set(gca,'YLim',[0 361+amp+0.2])
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
title([ '(c) Radial RFs for Station ' char(staname{1})],'FontSize',20 )

for j = 0:360
    time=[];
   drawingindex=(find(j==bazShow(:)));
   if (length(drawingindex)==0)
       continue;
   end
   mptm = find(RFs_Zam(drawingindex,:)==max(RFs_Zam(drawingindex,:)));
   time = (0:npts-1)/samprate - mptm(1)/samprate;
   tt= find( time >= tlength );   
   tmpZm=[]; tmpRm=[];
   tmpZm=RFs_Zam(drawingindex,:);
   tmpRm=RFs_Ram(drawingindex,:);         
   tmpZm=filtfilt(b,a,tmpZm').*Tf;
   tmpRm=filtfilt(b,a,tmpRm').*Tf;
   mptm1 = find(tmpZm(drawingindex,:)==max(tmpZm(drawingindex,:)));
   tmpRm=amp*tmpRm/max(abs(tmpRm));
   tmpZm=amp*tmpZm/max(abs(tmpZm));
    if (tmpRm(mptm1)>=0)
        ;
    else
          tmpRm=-tmpRm;
    end
    RFs_Rm_pick(j+1,:)=tmpRm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Drawing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   rr1m = tmpRm +j; rr1m(tt) = j;
   rr1m(1)=j; rr1m(end)=j;
   fii1m = find(rr1m <= j); frr1m = rr1m; frr1m(fii1m) = j;
   fii2m = find(rr1m >= j); frr2m = rr1m; frr2m(fii2m) = j;
   plot(time,rr1m,'k','LineWidth',1);   
   fill(time,frr1m','r','facealpha',0.7);
   fill(time,frr2m','b','facealpha',0.7);
% if Pick_tick~=1;
%    text(-3.0, j+0.2, [num2str(round(baz(baz_ind(j))))],'FontSize',15)
% else
%    text(-3.0, j+0.2, [num2str(bazShow(j))],'FontSize',15)
%    text(-2.2, j+0.2, num2str(j));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Picking the arrivals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if Pick_tick==1
   tt1=find(time>= PmsWindow(1) & time<= PmsWindow(2));
   Tpick(drawingindex)=find(tmpRm(tt1)==min(tmpRm(tt1)));
   plot(time(tt1(Tpick(drawingindex))),tmpRm(tt1(Tpick(drawingindex)))+j+0.01,'ko','linewidth',2);
   end
   drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stacking based on backazimuth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if QuadrantStack_tick==1;
if bazShow(drawingindex)<=90 & bazShow(drawingindex)>=0 %  Quadrant1 ?0-90?
   stackedRm_Q1=stackedRm_Q1+tmpRm;
   qua1_num=qua1_num+1;
elseif bazShow(drawingindex)<=180 & bazShow(drawingindex)>=90 %  Quadrant2 ?90-180?
   stackedRm_Q2=stackedRm_Q2+tmpRm;
   qua2_num=qua2_num+1;
elseif bazShow(drawingindex)<=270 & bazShow(drawingindex)>=180 %  Quadrant3 ?180-270?
   stackedRm_Q3=stackedRm_Q3+tmpRm;
   qua3_num=qua3_num+1;
elseif bazShow(drawingindex)<=360 & bazShow(drawingindex)>=270 %  Quadrant4 ?270-360?
   stackedRm_Q4=stackedRm_Q4+tmpRm;
   qua4_num=qua4_num+1;
end
end
stackedRm = stackedRm+tmpRm; % full stacked
end % end of loop for each event
stackedRm_Q1=stackedRm_Q1/qua1_num;
stackedRm_Q2=stackedRm_Q2/qua2_num;
stackedRm_Q3=stackedRm_Q3/qua3_num;
stackedRm_Q4=stackedRm_Q4/qua4_num;
stackedRm = stackedRm/length(bazShowind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Drawing Stacked RFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rr1m = stackedRm +j+1; rr1m(tt) = j+1;
rr1(1)=j+1; rr1(end)=j+1;
fii1 = find(rr1m <= j+1); frr1m = rr1m; frr1m(fii1) = j+1;
fii2 = find(rr1m >= j+1); frr2m = rr1m; frr2m(fii2) = j+1;
mptm = find(stackedRm==max(stackedRm));
time = (0:npts-1)/samprate - mptm/samprate;
plot(time, stackedRm+j+1,'b','LineWidth',1);
fill(time,frr1m','g');
fill(time,frr2m','w');
for tt=1:15
plot([tt tt],[0 361],'k--','LineWidth',1);
end
plot([0 0],[0 361],'k-','LineWidth',2);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Record the arrivals and the range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tt1=find(time>= PmsWindow(1) & time<= PmsWindow(2));
Tpick_stack=find(stackedRm(tt1)==min(stackedRm(tt1)));
for pick_i=1:length(Tpick)
    if (time(tt1(Tpick_stack))-uncertainty<=time(tt1(Tpick(pick_i))) & ...
        time(tt1(Tpick(pick_i)))<=time(tt1(Tpick_stack))+uncertainty)
       plot(time(tt1(Tpick(pick_i))),RFs_Rm_pick(pick_i,tt1(Tpick(pick_i)))+bazShow(pick_i)+0.01,'g*','linewidth',2);
    else
       Tpick(pick_i)=nan;
    end
end
plot(time(tt1(Tpick_stack)),stackedRm(tt1(Tpick_stack))+length(Tpick)+1+0.01,'ko','linewidth',2);
give_up=0;
Save_pick=input('\n Save the Pick (Yes,1/No,0)?');
if Save_pick==1
Save_Tpick=Tpick(~isnan(Tpick));
min_ar=0;max_ar=0;
min_ar=min(time(tt1(Save_Tpick)));
max_ar=max(time(tt1(Save_Tpick)));
else
give_up=input('\n give up?');
if give_up~=1
del_ra=input('\n Input dele #?');
res=setdiff(Tpick,Tpick(del_ra));
Save_Tpick=res(~isnan(res));
min_ar=0;max_ar=0;
min_ar=min(time(tt1(Save_Tpick)));
max_ar=max(time(tt1(Save_Tpick)));
end
end
if give_up~=1
fprintf(fid,'%s %2.2f %2.2f %2.2f %2.2f %2.2f\n',station,stalat(1),stalon(1),min_ar,max_ar,time(tt1(Tpick_stack)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % uncomment if you want to select
next=1 %input('\n Next Station (Yes,1/No,0)?')';
if next==1
    redraw=1;
%input('\n Save the stack (Yes,1/No,0)?')' % uncomment if you want to select
    Save= 0 %input('\n Save the stack (Yes,1/No,0)?')';
    if Save==1 & give_up~=1
    fname = ['RFs_stack' station '.mat']  
    maxbaz=max(baz); 
    minbaz=min(baz);
    eval(['save tmp_rf_moveout.mat stackedRm stackedRm_Q1 stackedRm_Q2 stackedRm_Q3 stackedRm_Q4 min_ar max_ar station nevents maxbaz minbaz stalat stalon'])
    unix(['mv tmp_rf_moveout.mat ' Resultfolder 'Stacked_Quadrant/' fname]);
    else
    end
    close all;
    continue;
end
close all
end

end
fclose(fid);

