

clear all; close all


datafolder = './RFs/';

stalist = textread([datafolder 'station_list_new'],'%s');

nstn = length(stalist)
% fid = fopen('PmS_MLD_manual.txt','w+');

for ss = 1:nstn, ss
station = char(stalist(ss))   

filename = dir( [datafolder '/' station '/RFs_waterlevel_' station '_*.mat'] );

maxpt = 4000;
  
RFs_Z=[]; RFs_R=[]; staname={}; baz=[]; iorid=[];stalat=[];stalon=[];
for i = 1:length(filename)
   files = load( [datafolder '/' station '/' filename(i).name] );

   RFs_Z(i,:) = files.RF_Zt(1:maxpt,:);
   RFs_R(i,:) = files.RF_Rt(1:maxpt,:);
   
   staname{i} = files.station; 
%    stalat(i) = files.stalat;
%    stalon(i) = files.stalon;
   baz(i) = files.backazi;
   iorid(i) = files.iorid;
   samprate = files.samprate;
end

dt=1/samprate;
npts=size(RFs_Z,2);   
nevents = length(baz)
[ixx,jj]  = sort(baz);

fband = [0.1 1.0];
Wn=fband.*(dt*2); % .1-1hz  , try different frequencies
[b,a]=butter(2,Wn);  % create filter coefficients (second order)
r=5; % secs
r=r/(dt*.5*npts);
Tf=tukeywin(npts,r);

figure('Position',[200 400 1000 1600]), hold on, box on, grid on

% subplot(1,2,1), hold on, box on, grid on
% xlabel('time after P-wave (seconds)', 'FontSize',16)
% set(gca,'xtick',[-10:2:20],'ygrid','on')
% set(gca,'XLim',[-5 10])
% set(gca,'YLim',[0 nevents+1])
% title([ 'Vertical RFs for Station ' char(staname{1})], 'FontSize',18)
% set(gca,'FontSize',14);

% subplot(1,2,2), hold on, box on, grid on
xlabel('time after P-wave (s)', 'FontSize',14)
set(gca,'Xtick',[-10:10:120],'ygrid','on')
set(gca,'XLim',[-1 10])
set(gca,'YLim',[0.5 nevents+2])
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
%text(0.4, nevents+1.5,'Stacked RF','FontSize',14,'Color','b');
title([ '(a) Radial RFs for Station ' char(staname{1})],'FontSize',18 )
set(gca,'FontSize',12);

stackedR = zeros(maxpt,1);
for j = 1:nevents
    
   mpt = find(RFs_Z(jj(j),:)==max(RFs_Z(jj(j),:)));
   if mpt~=601, disp('mpt != 601'), end
   
   time = (0:npts-1)/samprate - mpt/samprate;
   tt= find( time >= 115 );
   
   tmpZ=[]; tmpR=[];
   
   tmpZ=RFs_Z(jj(j),:);
   tmpR=RFs_R(jj(j),:);
   
%    if RFs_R(jj(j),mpt)>=0, 
%    tmpR=RFs_R(jj(j),:);
%    elseif RFs_R(jj(j),mpt)<0, 
%    tmpR=-RFs_R(jj(j),:); 
%    end
      
   tmpZ=filtfilt(b,a,tmpZ').*Tf;
   tmpR=filtfilt(b,a,tmpR').*Tf;
   
   tmpR=tmpR/max(tmpR);
   tmpZ=tmpZ/max(tmpZ);
   
   stackedR = stackedR+tmpR;
   
   %%%%% vertical RFs
%    rr1 = tmpZ +j; rr1(tt) = j;
%    rr1(1)=j; rr1(end)=j;
%    fii1 = find(rr1 <= j); frr1 = rr1; frr1(fii1) = j;
%    fii2 = find(rr1 >= j); frr2 = rr1; frr2(fii2) = j;
   
%    subplot(1,2,1)
%    plot(time,rr1,'k'); 
%    fill(time,frr1',[0.5 0.5 0.5]);
%    fill(time,frr2','w')  
%    
%    text(-4.8, j+0.3, num2str(iorid(jj(j))) );
%    text(8.0, j+0.3, [num2str(round(baz(jj(j))))])
   
   %%%%% Radial RFs
   rr1 = tmpR +j; rr1(tt) = j;
   rr1(1)=j; rr1(end)=j;
   fii1 = find(rr1 <= j); frr1 = rr1; frr1(fii1) = j;
   fii2 = find(rr1 >= j); frr2 = rr1; frr2(fii2) = j;
   
%    subplot(1,2,2)
   plot(time,rr1,'b'); 
   fill(time,frr1','r');
   fill(time,frr2','g')  

%    text(-2.2, j+0.3, num2str( iorid(jj(j)) ));
   text(-1.5, j+0.2, [num2str(round(baz(jj(j))))],'FontSize',16)
    
   drawnow
end % end of loop for each event
stackedR = stackedR/nevents;

rr1 = stackedR +nevents+1; rr1(tt) = nevents+1;
rr1(1)=nevents+1; rr1(end)=nevents+1;
fii1 = find(rr1 <= nevents+1); frr1 = rr1; frr1(fii1) = nevents+1;
fii2 = find(rr1 >= nevents+1); frr2 = rr1; frr2(fii2) = nevents+1;
   
%    subplot(1,2,2)
plot(time,rr1,'b'); 
fill(time,frr1','b');
fill(time,frr2','w')  

plot(time, stackedR+nevents+1,'b','LineWidth',2);

for tt=1:120
plot([tt tt],[0 nevents+1],'k--','LineWidth',1);
end
% for tt=0:5:20
% plot([tt tt],[0 nevents+1],'k-','LineWidth',2);
% end
plot([0 0],[0 nevents+1],'k-','LineWidth',2);

% set(gcf, 'PaperPositionMode','auto');
% figname = ['RFs_' char(station) '_' num2str(fband(1)) '-' num2str(fband(2)) ];
% eval(['print -depsc ' figname '.eps']);
disp('Enter for next station ...')
pause
close all;

%%% manually pick up PmS, and MLD
% [pms, amp] = ginput;
% 
% fprintf(fid,'%f %f %s %f %f \n',stalon(j), stalat(j), char(station), pms(1), pms(2)); 
% disp([j ' Click Enter for next event ...']);

end
% fclose(fid);

