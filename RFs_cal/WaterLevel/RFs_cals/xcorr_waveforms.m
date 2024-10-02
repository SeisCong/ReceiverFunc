%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% THIS CODE IS USED FOR INTERACTIVE CROSS-CORRELATION AND PLOT FOR EACH EVENT
%
%
% Original one from Haiying Gao, modified by Cong Li
% 
% 19/5/16 Add the direction correction for the eastern, northern components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

samprate=20;                % desired sample rate for data
count_c=1;istnn_c={};
fid1 = fopen('AngCorrection.txt', 'w');
fid2 = fopen('AngCorrection_Station.txt', 'w');
fid3 = fopen('AngCorrection_delete.txt', 'w');
filename = dir('./data/M5.0/*.mat');
    
for fi =1:length(filename), fi, filename(fi).name
    
    HZ=[]; H1=[]; H2=[]; istn={}; baz=[]; dep=[]; slowP=[];
    MagB=[]; delta=[]; ymdtime={}; slat=[]; slon=[]; olat=[]; olon=[];
    hangH1=[];hangH2=[];E={};N={}; Rt={};Tt={}; record=[];Zc=[];Rc=[];Tc=[];
    Z1=[];R1=[];T1=[];
    files = load( ['./data/M5.0/' filename(fi).name] ); 
    
    eorid = filename(fi).name(5:end-4)
    
    Z = files.Zcomp;
    H1r = files.H1comp;
    H2r = files.H2comp;
    iorid = files.iorid;
    istn = files.istn;
    baz = files.baz;
    dep = files.dep;
    slowP = files.slowP;
    MagB = files.MagB;
    delta = files.delta;
    ymdtime = files.ymdtime;
    slat = files.slat;
    slon = files.slon;
    olat = files.olat;
    olon= files.olon;
    hangH1=files.Ehang;
    hangH2=files.Nhang;
    sps=files.sampleo;
    
    nStn = length(istn);     % number of stations for this event
    nt = size(Z,1);
    Zt = Z;
    % Correction of direction for those non-eastern and non-northern record
    plottick=1;
    count_plot=1;
        for iar=1:nStn
            Er=[];Nr=[];hangE=0;hangN=0;
        if (abs(hangH1(iar)-90)<=45||abs(hangH1(iar)-270)<=45) % Judging the near east-west and north-south component
            Er=H1r{1,iar};
            Nr=H2r{1,iar};
            hangE=hangH1(iar);
            hangN=hangH2(iar);
        else
            Er=H2r{1,iar};
            Nr=H1r{1,iar};
            hangE=hangH2(iar);
            hangN=hangH1(iar);
        end
         if (abs(hangE-90)<=0.01&&(abs(hangN-360)<=0.01||abs(hangN<=0.01))) % Normal, tolerant error is equal to 0.01
          E{1,iar}=Er;
          N{1,iar}=Nr;
         elseif (abs(hangE-hangN)==90 || abs(hangE-hangN)==270)      
          E{1,iar}=Er*sind(hangE)+Nr*sind(hangN);
          N{1,iar}=Er*cosd(hangE)+Nr*cosd(hangN);
          record(count_plot)=iar;
          istnn_c{count_c}=sprintf('%s  %2.2f  %2.2f ',istn{iar},hangE,hangN);
          fprintf(fid1,'%s  %2.2f  %2.2f %d %s\n',istn{iar},hangE,hangN,fi,filename(fi).name); 
          count_c=count_c+1;
          count_plot=count_plot+1;
          plottick=1;
         else
          E{1,iar}=zeros(130*sps(iar),1);
          N{1,iar}=zeros(130*sps(iar),1);
          fprintf(fid3,'%s  %2.2f  %2.2f %d %s\n',istn{iar},hangE,hangN,fi,filename(fi).name);
        end
        end
    %%rotate ENZ to RTZ components
    for iar=1:nStn    
      ang=270-baz(iar);
      Rt{1,iar} = E{1,iar}*cosd(ang)+N{1,iar}*sind(ang);
      Tt{1,iar} = -E{1,iar}*sind(ang)+N{1,iar}*cosd(ang);
    end

% =========================================================================
%%      Plot for the corrected and rotated component
% =========================================================================
% if (plottick==1)
%   close all;
%   figure('Position',[400 400 1400 1200]);
%   offset=0;
%   subplot(1,2,1);
%   for ind=record(1:end) 
%      if (ind==record(1))
%          ;
%      else
%          offset=offset+2;
%      end
%      plot(H1r{1,ind}/max(H1r{1,ind})+offset,'r')
%      offset=offset+2;
%      hold on;
%      plot(E{1,ind}/max(E{1,ind})+offset,'b')
%      offset=offset+2;
%      hold on;
%      plot(Rt{1,ind}/max(Rt{1,ind})+offset,'k')
%      text(130*sps{ind},offset,istn{ind},'Fontsize',10)  
%   end
%   legend({'Raw east component','Corrected east component','R component'},'FontSize',15,'Location','southeast')
%   title(['Correction and rotation of East component for event: ' num2str(eorid)],'FontSize',15 )  
%   subplot(1,2,2);
%   offset=0;
%   for ind=record(1:end)
%      if (ind==record(1))
%          ;
%      else
%          offset=offset+2;
%      end
%      plot(H2r{1,ind}/max(H2r{1,ind})+offset,'r')
%      offset=offset+2;
%      hold on;
%      plot(N{1,ind}/max(N{1,ind})+offset,'b')
%      offset=offset+2;
%      hold on;
%      plot(Tt{1,ind}/max(Tt{1,ind})+offset,'k')
%      text(130*sps{ind},offset,istn{ind},'Fontsize',10)  
%   end
%   legend({'Raw north component','Corrected north component','T component'},'FontSize',15,'Location','southeast')
%   title(['Correction and rotation of North component for event: ' num2str(eorid)],'FontSize',15 )
% end
% =========================================================================
   
   Za=Zt;
   Ra=Rt;
   Ta=Tt;
   
   close all;
   figure('Position',[400 400 1400 1200]);
      
   H1=subplot('position',[0.1 0.03 0.6 0.9]); cla, hold on   % traces
   H2=subplot('position',[0.75 0.03 0.15 0.9]); cla, hold on   % delay dots

   lagTotal=zeros(nStn,1);
   StnFlag=ones(nStn,1);   % 1 is good, 0 is out
   selv_corr=zeros(nStn,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Plot the raw vertical component for all stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    subplot(H1); 
%    offset=0;
%    for ind=1:nStn 
%      plot(Za{1,ind}/max(Za{1,ind})+offset)
%      offset=offset+1;
%    end
%    title( ['vertical motion for event: ' num2str(eorid)] )
%    axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Zb=Za;
   Zcc=Zb;
  
   Rb=Ra;
   Rcc=Rb;
   Tb=Ta;
   Tcc=Tb;

   dt=1/samprate;
   Wn=[.05,1.2].*(dt*2); % .1-1hz  , try different frequencies
   [b,a]=butter(2,Wn);  % create filter coefficients (second order)
   r=8; % secs 

   for ff=1:nStn
   %% Resampling
   Z1(:,ff)=resample(Zb{1,ff},samprate, sps(ff)); % resampled raw data 
   R1(:,ff)=resample(Rb{1,ff},samprate, sps(ff));
   T1(:,ff)=resample(Tb{1,ff},samprate, sps(ff));
   Zcc{1,ff}=resample(Zcc{1,ff},samprate, sps(ff)); 
   Rcc{1,ff}=resample(Rcc{1,ff},samprate, sps(ff));
   Tcc{1,ff}=resample(Tcc{1,ff},samprate, sps(ff));  
   
   ntc=size(Zcc{1,ff},1);  
   r=r/(dt*.5*ntc);
   Tf=tukeywin(ntc,r);
   Zcc{1,ff}=filtfilt(b,a,Zcc{1,ff}.*Tf); % resampled filtered data
   Rcc{1,ff}=filtfilt(b,a,Rcc{1,ff}.*Tf);
   Tcc{1,ff}=filtfilt(b,a,Tcc{1,ff}.*Tf);   
  
   Zc(:,ff)=Zcc{1,ff};
   Rc(:,ff)=Rcc{1,ff};
   Tc(:,ff)=Tcc{1,ff};
   end 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%       PLOT FILTERED 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    close all;     
%    figure('Position',[400 400 800 1200]); 
%    offset=0;
%    for ii=1:nStn
%     mpt = find(Zcc{1,ii}==max(Zcc{1,ii}));
%      cc1=mpt-29*samprate;
%     time=(1:length(Zcc{1,ii}(cc1:end)))/samprate-29;
% %     time=(1:length(Zcc{1,ii}))/samprate-30;
%     subplot(1,2,1);%,'position',[0.1 0.1 0.5 0.8]);    
%     plot(time,Zcc{1,ii}(cc1:end)/max(Zcc{1,ii}(cc1:end))+offset,'k','LineWidth',1);
% %     plot(time,Zcc{1,ii}/max(Zcc{1,ii})+offset,'k','LineWidth',1);
%     hold on;
%     xlim([-2,10]);
%     ylim([-2,55]);
%     set(gca,'fontsize',18);
%     xlabel('Time after P arrival (s)','Fontsize',20);
%     text(-2.2,offset,istn{ii},'Fontsize',15,'FontWeight','bold','HorizontalAlignment','right')
%     %text(0,-5.0,'Filtered Z component','Fontsize',20,'FontWeight','bold')
%     title('(a) Vertical Component','Fontsize',30)
%     box on, grid on
%     ax=gca;
%     ax.XTick=[-2:2:16]
%     ax.YTick=[-10:10:-10]
%     ax.GridAlpha = 0.5;
%     subplot(1,2,2);
%     plot(time,Rcc{1,ii}(cc1:end)/max(Rcc{1,ii}(cc1:end))+offset,'k','LineWidth',1);
% %     plot(time,Rcc{1,ii}/max(Rcc{1,ii})+offset,'k','LineWidth',1);
%     hold on;
%     xlim([-2,10]);
%     ylim([-2,55]);
%     set(gca,'fontsize',18);
%     xlabel('Time after P arrival (s) ','Fontsize',20);
%     text(-2.2,offset,istn{ii},'Fontsize',15,'FontWeight','bold','HorizontalAlignment','right')
%     %text(0,-5.0,'Filtered R component','Fontsize',20,'FontWeight','bold')
%     title('(b) Radial Component','Fontsize',30)
%     box on, grid on
%     ax=gca;
%     ax.XTick=[-2:2:16]
%     ax.YTick=[-10:10:-10]
%     offset=offset+2;
%     ax.GridAlpha = 0.5;
%    end
% %     title('Filtered Waveform');
%     pause;
%     set(gcf,'PaperPositionMode','auto')
%     figname = ['VandRComponent'];
%     eval(['print -depsc ' figname '.eps']);
%     close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the filtered vertical component 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         H1 = resample(H1,samprate, sps);
%         H2 = resample(H2,samprate, sps);
%    figure('Position',[400 400 1400 1200]);
%    offset=0;
%    for ind=1:nStn 
%      plot(Zc(:,ind)/max(Zc(:,ind))+offset);
%      hold on;
%      offset=offset+2;
%    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   pre=1;
   post=ntc;
   tx=[];

   Tcol=0;                 % highlighted trace
   extraT=samprate*5;


for kount=1:1000 % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
 if kount==1 %------------------------------------- %
   hh=6;                                            %
%  elseif kount==2                                    %
%    hh=8;                                            %
 else                                               %
   hh=input(['\n\n-Your choice:' ...                %
   '\n  (1) add (or remove) time line'...           %
   '\n  (2) filter traces' ...                      %
   '\n  (3) manually time shift a trace' ...        %
   '\n  (4) cross correlate traces' ...             %
   '\n  (5) flip a trace' ...                       %
   '\n  (6) reset all lag times to zero' ...        %
   '\n  (7) eliminate a trace' ...                  %
   '\n  (8) pick phase to corect arrival times' ... %
   '\n  (9) highlight a trace' ...                  %
   '\n  (0) continue to next event' ...             %
   '\n==> ']);                                      %
 end %----------------------------------------------%
 
 if hh==0 % <<<<<<<<<<<<<<< go to next event >>>>>>>>>>>>>>> (0)
   break;
 elseif hh==1 % <<<<<<<<<<< add time line <<<<<<<<<<<<<<<<<< (1)
   t1=input('time (neg to clear): ');
   if t1>0
     tx=[tx t1];
   else
     tx=[];
   end
 elseif hh==2 % <<<<<<<<<<< re-filter traces <<<<<<<<<<<<<<< (2)
     cf=input('-center frequency (Hz): ');
     if isempty(cf), cf=1; end
     hw=input('-half-width (Hz): ');
     if isempty(hw), hw=0.5; end
     A = Gfilt(cf,hw,ntc,samprate);
 %    Zc=Zb;
     for ff=1:nStn
     Zc(:,ff)=ifft(fft(Zc(:,ff)).*A);
     Rc(:,ff)=ifft(fft(Rc(:,ff)).*A);
     Tc(:,ff)=ifft(fft(Tc(:,ff)).*A);
     end      

 elseif hh==3  % <<<<<<<<<<< manually shift traces <<<<<<<<< (3)
   Tnum=input('-enter trace number to be shifted: ');
     if Tnum<1 || Tnum>nStn
       Tnum=input('-BAD TRACE NUMBER!  Enter trace number to be shifted: ');
     end
   Tshift=input('\n-enter shift (in points): ');
   lagTotal(Tnum) = lagTotal(Tnum) + Tshift;
%    lagAvg   = sum(StnFlag.*lagTotal)/sum(StnFlag);
%    lagTotal = lagTotal-lagAvg;
 elseif hh==4  % <<<<<<<<<<< cross correlate <<<<<<<<<<<<<<< (4)
   pre1=pre;  %%%% add by cong
   post1=post; %%%% add by cong
   pre =input('\n-Enter your desired X-corr start point: ');
   post=input('\n-Enter your desired X-corr end point: ');
   laglim = input('\n-Enter lag limit (in points): ');
   if isempty(pre),    pre =pre1;  end % Default pre==1
   if isempty(post),   post=post1; end % Default post==end
   if isempty(laglim), laglim=0.75*samprate; end
   pre1=pre; post1=post;
   if laglim<0
     laglim = -laglim;
     extraT = input('\n-Enter plot time lead (in points): ');
   end
   Zd=zeros(ntc,nStn);  % Zd is shifted traces for X-corr
   for k=1:nStn
     lagT=round(lagTotal(k));
     if lagT<0
       Zd(1:ntc+lagT,k)=Zc(1-lagT:ntc,k);
       Rd(1:ntc+lagT,k)=Rc(1-lagT:ntc,k);
       Td(1:ntc+lagT,k)=Tc(1-lagT:ntc,k);            
     else
       Zd(1+lagT:ntc,k)=Zc(1:ntc-lagT,k);
       Rd(1+lagT:ntc,k)=Rc(1:ntc-lagT,k);
       Td(1+lagT:ntc,k)=Tc(1:ntc-lagT,k);
     end
   end
   lagtime  =  LSQxcorr(Zd,pre,post,laglim);
   lagTotal = lagTotal+lagtime;
   lagAvg   = sum(StnFlag.*lagTotal)/sum(StnFlag);
   lagTotal = lagTotal-lagAvg;
 elseif hh==5 % <<<<<<<<<<< flip trace <<<<<<<<<<<<<< (5)
     Tnum=input('-enter trace number to be flipped: ');
     if Tnum<1 || Tnum>nStn
       Tnum=input('-BAD TRACE NUMBER!  Enter trace number to be flipped: ');
     end
     Zc(:,Tnum) = -Zc(:,Tnum);
     Zb(:,Tnum) = -Zb(:,Tnum);
     Za(:,Tnum) = -Za(:,Tnum);
     Rc(:,Tnum) = -Rc(:,Tnum);
     Rb(:,Tnum) = -Rb(:,Tnum);
     Ra(:,Tnum) = -Ra(:,Tnum);
     Tc(:,Tnum) = -Tc(:,Tnum);
     Tb(:,Tnum) = -Tb(:,Tnum);
     Ta(:,Tnum) = -Ta(:,Tnum);         
 elseif hh==6 % <<<<<<<<<<< reset lag times to zero <<<<<<<< (6)
   lagTotal = 0*lagTotal;
 elseif hh==7 % <<<<<<<<<<< remove trace <<<<<<<<<<<<<<<<<<< (7)
   StnOut=input('\n-Station number to be eliminated: ');
   StnFlag(StnOut)=0;
   lagTotal(StnOut)=0;
 
 elseif hh==9 % <<<<<<<<<<< highlight a trace <<<<<<<<<<<<<< (9)3
   Tcol=input('-enter trace number to be highlighted: ');
 end
 % _____ main plotting _______
   subplot(H1); cla, hold on, grid on
   bb = round(2+max(abs(lagTotal)));
   dd=round(delta*10)/10;
   bazAvg=mean(baz);
   for k=1:nStn  
     if StnFlag(k)==1
       col='b';
       if k==Tcol, col='g'; end
       traceZ = Zc(1:ntc,k);
       traceR = Rc(1:ntc,k);
       traceT = Tc(1:ntc,k);
       ss = round(lagTotal(k));
       Zx(bb:ntc-bb,k) = traceZ(bb-ss:ntc-bb-ss);   % bb is the standard time, ss is the relative time    
       Rx(bb:ntc-bb,k) = traceR(bb-ss:ntc-bb-ss);      
       Tx(bb:ntc-bb,k) = traceT(bb-ss:ntc-bb-ss);
       
       % unfiltered data (raw waveform)
        traceZr = Z1(:,k);
        traceRr = R1(:,k);
        traceTr = T1(:,k);
        ss = round(lagTotal(k));
        
        Zraw(:,k) = traceZr(bb-ss:ntc-bb-ss);
        Rraw(:,k) = traceRr(bb-ss:ntc-bb-ss);
        Traw(:,k) = traceTr(bb-ss:ntc-bb-ss); 
%        Zraw=Zb;
%        Rraw=Rb;
%        Traw=Tb;
 
       
       p1=max([    bb  pre-extraT]);
       if (hh==3)
       p2=min([ntc post+extraT]);
       else
       p2=min([ntc-bb post+extraT]);
       end
       traceZ = traceZ/(2*max(traceZ(p1-ss:p2-ss)));
       traceR = traceR/(2*max(traceR(p1-ss:p2-ss)));
       traceT = traceT/(2*max(traceT(p1-ss:p2-ss)));
       
       plot( p1:p2, traceZ(p1-ss:p2-ss)+k,col )
       text(p1+(p2-p1)*0.010,k+0.27, num2str(k) );
       text(p1+(p2-p1)*0.045,k+0.27, istn(k) );
       text(p1-0.08*(p2-p1),k, num2str(dd(k)) );
     else, continue
     end
   end
   plot([pre pre], [0 nStn+1],'k')
   plot([post post], [0 nStn+1],'k')
   set(H1,'YTickLabel',[])
   if isempty(tx)==0
     for kk=1:length(tx), plot([1 1]*tx(kk), [-1 nStn+1],'r'), end
   end

   ee = num2str(iorid(k));
   bb = num2str(round(bazAvg));
   title(['aligned  Event number ' ee '.   Baz= ' bb])
   axis tight
   ylim([-0.2 nStn+1])

 subplot(H2); cla, hold on, grid on
   plot([0 0], [0 nStn+1],'k')
   for k=1:nStn
     if StnFlag(k)==1
       ss = round(-lagTotal(k))-selv_corr(k);
       plot(ss,k,'r*','MarkerSize',10)
     end
   end
   xlim([min(-lagTotal)-samprate/5 max(-lagTotal)+samprate/5])
   ylim([-0.2 nStn+1])
   title('Delay Times')
   set(H2,'YAxisLocation','right')
   
end
   
% save data for receiver functions
 sdf=input(['\n\n-Your choice:' ...                %
    '\n  (0) save result for receiver functions' ...
    '\n  (1) unsave data and go to next event' ...
    '\n==> ']);
 if sdf ==0
     idreal = find(StnFlag~=0);
   if ~isempty(idreal)
      Zx_raw = Zraw(:,idreal); Rx_raw = Rraw(:,idreal); Tx_raw = Traw(:,idreal); %Traw(:,idreal);
     istn_new=istn(idreal); iorid_new=iorid(idreal); slowP_new=slowP(idreal); delta_new=delta(idreal);
     baz_new=baz(idreal); ymdtime_new=ymdtime(idreal); dep_new=dep(idreal);
     MagB_new=MagB(idreal); olat_new=olat(idreal); olon_new=olon(idreal);
     slat_new=slat(idreal); slon_new=slon(idreal); 
     
%      files = ['ZRT_' num2str(eorid) '.mat'];
%      save  ZRT_data.mat   Zx_raw Rx_raw Tx_raw  istn_new iorid_new slowP_new baz_new  delta_new ...
%         ymdtime_new dep_new MagB_new olat_new olon_new slat_new slon_new samprate   % save data for receiver functions
% 
%      unix(['mv ZRT_data.mat data/data_xcorr/'  files]);  
   end
       
 end %
 
end % loop for each event =================================================
AA=unique(istnn_c);
 for ii=1:size(AA,2)
 fprintf(fid2, '%s\n', AA{ii});
 end
fclose(fid1);
fclose(fid2);
fclose(fid3);
%=========================================================================
