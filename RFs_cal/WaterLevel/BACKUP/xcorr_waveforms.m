% Interactive cross-correlation and plot for each event.

clear all; close all;

samprate=20;                % desired sample rate for data
    
filename = dir('./data/M5.0/*.mat');
    
for fi = 49:length(filename), fi, filename(fi).name
    
    HZ=[]; H1=[]; H2=[]; istn={}; baz=[]; dep=[]; slowP=[];
    MagB=[]; delta=[]; ymdtime={}; slat=[]; slon=[]; olat=[]; olon=[];

    files = load( ['./data/M5.0/' filename(fi).name] ); 
    
    eorid = filename(fi).name(5:end-4)
    
    Z = files.Zcomp;
    Er = files.H1comp;
    Nr = files.H2comp;
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
    hangE=files.Ehang;
    hangN=files.Nhang;
    
    nStn = length(istn);     % number of stations for this event
    nt = size(Z,1);
    
    
    Zt = zeros(nt, nStn);
    Rt = zeros(nt, nStn);
    Tt = zeros(nt, nStn);
    E = zeros(nt, nStn);
    N = zeros(nt, nStn);
    Zt = Z;
    % Correction of direction for those non-eastern and non-northern record
    for iar=1:nStn
        if (hangE(iar)==90)
            continue;
        else
        Eang_corr=hangE(iar)-90; % clockwise is positive
        Nang_corr=hangN(iar)-360;
        E(:,iar)=Er(:,iar)*cosd(Eang_corr)+Nr(:,iar)*sind(Eang_corr);
        N(:,iar)=Nr(:,iar)*cosd(Eang_corr)-Er(:,iar)*sind(Eang_corr);
        if (hangN(iar)==90)
            N(:,iar)=-N(:,iar);
        end
        end
    end
        
    %%rotate ENZ to RTZ components
    for iar=1:nStn    
      ang=270-baz(iar);
      Rt(:,iar) = E(:,iar)*cosd(ang)+N(:,iar)*sind(ang);
      Tt(:,iar) = -E(:,iar)*sind(ang)+N(:,iar)*cosd(ang);
    end
        
% =========================================================================

   Za=zeros(nt,nStn);        % the raw Z seismograms for this event
   Ra=zeros(nt,nStn);
   Ta=zeros(nt,nStn); 
   
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

   subplot(H1); 
   offset=0;
   for ind=1:nStn 
     plot(Za(:,ind)/max(Za(:,ind))+offset)
     offset=offset+1;
   end
   title( ['vertical motion for event: ' num2str(eorid)] )
   axis tight

   Zb=Za;
   Zc=Zb;
   ntc=size(Zc,1);
   Rb=Ra;
   Rc=Rb;
   Tb=Ta;
   Tc=Tb;

   dt=1/samprate;
   Wn=[.3,1.5].*(dt*2); % .1-1hz  , try different frequencies
   [b,a]=butter(2,Wn);  % create filter coefficients (second order)
   r=8; % secs
   r=r/(dt*.5*ntc);
   Tf=tukeywin(ntc,r);
   for ff=1:nStn
   Zc(:,ff)=filtfilt(b,a,Zc(:,ff).*Tf);
   Rc(:,ff)=filtfilt(b,a,Rc(:,ff).*Tf);
   Tc(:,ff)=filtfilt(b,a,Tc(:,ff).*Tf);
   end
 
   pre=1;
   post=ntc;
   tx=[];

   Tcol=0;                 % highlighted trace
   extraT=samprate*5;
end
for kount=1:1000 % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
 if kount==1 %------------------------------------- %
   hh=6;                                            %
 elseif kount==2                                    %
   hh=8;                                            %
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
     Zc=Zb;
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
   lagAvg   = sum(StnFlag.*lagTotal)/sum(StnFlag);
   lagTotal = lagTotal-lagAvg;
 elseif hh==4  % <<<<<<<<<<< cross correlate <<<<<<<<<<<<<<< (4)
   pre =input('\n-Enter your desired X-corr start point: ');
   post=input('\n-Enter your desired X-corr end point: ');
   laglim = input('\n-Enter lag limit (in points): ');
   if isempty(pre),    pre =pre1;  end
   if isempty(post),   post=post1; end
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
 
 elseif hh==9 % <<<<<<<<<<< highlight a trace <<<<<<<<<<<<<< (9)
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
       Zx(bb:ntc-bb,k) = traceZ(bb-ss:ntc-bb-ss);       
       Rx(bb:ntc-bb,k) = traceR(bb-ss:ntc-bb-ss);      
       Tx(bb:ntc-bb,k) = traceT(bb-ss:ntc-bb-ss);
       
       % unfiltered data (raw waveform)
       traceZr = Zb(1:ntc,k);
       traceRr = Rb(1:ntc,k);
       traceTr = Tb(1:ntc,k);
       ss = round(lagTotal(k));
       Zraw(bb:ntc-bb,k) = traceZr(bb-ss:ntc-bb-ss);
       Rraw(bb:ntc-bb,k) = traceRr(bb-ss:ntc-bb-ss);
       Traw(bb:ntc-bb,k) = traceTr(bb-ss:ntc-bb-ss);  
       
       p1=max([    bb  pre-extraT]);
       p2=min([ntc-bb post+extraT]);
       traceZ = traceZ/max(traceZ(p1-ss:p2-ss));
       traceR = traceR/max(traceR(p1-ss:p2-ss));
       traceT = traceT/max(traceT(p1-ss:p2-ss));
       
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
       plot(ss,k,'r.')
     end
   end
   xlim([min(-lagTotal)-samprate/5 max(-lagTotal)+samprate/5])
   ylim([-0.2 nStn+1])
   title('Delay Times')
   set(H2,'YAxisLocation','right')
   
end
   
% save data for receiver functions
% %  sdf=input(['\n\n-Your choice:' ...                %
% %     '\n  (0) save result for receiver functions' ...
% %     '\n  (1) unsave data and go to next event' ...
% %     '\n==> ']);
% %  if sdf ==0
% %      idreal = find(StnFlag~=0);
% %    if ~isempty(idreal)
% %      Zx_raw = Zraw(:,idreal); Rx_raw = Rraw(:,idreal); Tx_raw = Traw(:,idreal);
% %      istn_new=istn(idreal); iorid_new=iorid(idreal); slowP_new=slowP(idreal); delta_new=delta(idreal);
% %      baz_new=baz(idreal); ymdtime_new=ymdtime(idreal); dep_new=dep(idreal);
% %      MagB_new=MagB(idreal); olat_new=olat(idreal); olon_new=olon(idreal);
% %      slat_new=slat(idreal); slon_new=slon(idreal); 
% %      
% %      files = ['ZRT_' num2str(eorid) '.mat'];
% %      save  ZRT_data.mat   Zx_raw Rx_raw Tx_raw  istn_new iorid_new slowP_new baz_new  delta_new ...
% %         ymdtime_new dep_new MagB_new olat_new olon_new slat_new slon_new samprate   % save data for receiver functions
% % 
% %      unix(['mv ZRT_data.mat data_xcorr/M5.4_7.0/'  files]);  
% %    end
% %         
% %  end %
% %  
% % end % loop for each event =================================================
% =========================================================================
