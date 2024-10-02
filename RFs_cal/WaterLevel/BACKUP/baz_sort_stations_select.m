

clear all; close all


datafolder = './RFs/';

stalist = textread([datafolder 'station_list'],'%s');

nstn = length(stalist)


maxpt = 4000;
for ss = 1:nstn, ss
station = char(stalist(ss))   

filename = dir( [datafolder '/' station '/RFs_waterlevel_' station '_*.mat'] );
  
RFs_Z=[]; RFs_R=[]; staname={}; baz=[]; iorid=[];stalat=[];stalon=[];
for i = 1:length(filename)
   files = load( [datafolder '/' station '/' filename(i).name] );

   RFs_Z(i,:) = files.RF_Zt(1:maxpt,:);
   RFs_R(i,:) = files.RF_Rt(1:maxpt,:);
   
   staname{i} = files.station; 

   baz(i) = files.backazi;
   iorid(i) = files.iorid;
   samprate = files.samprate;
end
dt=1/samprate;
npts=size(RFs_Z,2);   
nevents = length(baz);
[ixx,jj]  = sort(baz);

fband = [0.2 1.2];
Wn=fband.*(dt*2); % .1-1hz  , try different frequencies
[b,a]=butter(2,Wn);  % create filter coefficients (second order)
r=5; % secs
r=r/(dt*.5*npts);
Tf=tukeywin(npts,r);

n =50;
for cc=1:n:nevents

figure('Position',[200 400 1000 1600]), hold on, box on, grid on

xlabel('time after P-wave (s)', 'FontSize',14)
set(gca,'Xtick',[-10:5:10],'ygrid','on')
set(gca,'XLim',[-1 10])
set(gca,'YLim',[cc-1 cc+n])
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
title([ '(a) Radial RFs for Station ' char(staname{1})],'FontSize',18 )
set(gca,'FontSize',12);

for j = cc:cc+n-1    
    
   mpt = find(RFs_Z(jj(j),:)==max(RFs_Z(jj(j),:)));
   if mpt~=601, disp('mpt != 601'), end
      
   time = (0:npts-1)/samprate - mpt/samprate;
   tt= find( time >= 15 );
   
   tmpZ=[]; tmpR=[];
   
   tmpZ=RFs_Z(jj(j),:);
%     tmpR=RFs_R(jj(j),:);
      
   if RFs_R(jj(j),mpt)>0, 
   tmpR=RFs_R(jj(j),:);
   elseif RFs_R(jj(j),mpt)<=0, 
   tmpR=-RFs_R(jj(j),:); 
   end
         
   tmpZ=filtfilt(b,a,tmpZ').*Tf;
   tmpR=filtfilt(b,a,tmpR').*Tf;
   
   tmpR=tmpR/max(tmpR);
   tmpZ=tmpZ/max(tmpZ);
         
   %%%%% Radial RFs
   rr1 = tmpR +j; rr1(tt) = j;
   rr1(1)=j; rr1(end)=j;
   fii1 = find(rr1 <= j); frr1 = rr1; frr1(fii1) = j;
   fii2 = find(rr1 >= j); frr2 = rr1; frr2(fii2) = j;
   
%    subplot(1,2,2)
   plot(time,rr1,'b'); 
   fill(time,frr1','r');
   fill(time,frr2','g')  

   text(-1.5, j+0.3, num2str( iorid(jj(j)) ));
   text(-2.5, j+0.2, [num2str(round(baz(jj(j))))],'FontSize',16)
    
   drawnow
end % end of loop for each event

for tt=1:10
plot([tt tt],[0 nevents+1],'k--','LineWidth',1);
end

plot([0 0],[0 nevents+1],'k-','LineWidth',2);

pause
close all
end

end

