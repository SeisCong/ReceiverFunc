% this program uses least square method to estimate the optimitized
% relative arrival time, which will give the lag times
function [lagtime] = LSQxcorr(Ev,t1,t2,laglim)
% to time shifts applied
% [t1,t2] interval used for cross correlation
% laglim max lag

nStn = size(Ev,2);
xc = zeros(nStn,nStn);
nt=t2-t1+1;
i2nt=1:nt;
Ex=zeros(nt+2*laglim,nStn);
for k =1:nStn
Ex(laglim+1:nt+laglim,k) = Ev(t1:t2,k);
end
for i=1:nStn-1
  shifted=Ex(:,i);
  for j=i+1:nStn
      fixed = Ex(:,j);
      for lag=-laglim:laglim
       aa(lag+laglim+1) = sum(fixed(laglim+i2nt).*shifted((laglim-lag)+i2nt));
      end
      [mm,xc(i,j)]=max(aa);
  end 
end 

for k1=1:nStn-1
    for k2=k1+1:nStn
    xc(k1,k2) = xc(k1,k2)-laglim-1;
    end
end

xc = xc-xc';

ave = sum(xc,2)/nStn;
lagtime = ave;

