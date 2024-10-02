function [rftn] = dcv(Z,R,sampling_rate,Ns, Ne,nTap,NW,time_shift)

% function [rftn] = dcv(Z,R,sampling_rate,water_level,time_shift) 
%
% This is a script for receveiver function deconvolution.
% via the water level method.
% 
% -Z is the Z component.
% -R is either the radial or transverse component.
% -sampling_rate is the sampling rate.
% -water_level is the deconvolution parameter of the same name.
% -time_shift is the amount of time in the signal before the P arrival.
%


original_length = length(R);
fft_points = 2^(ceil(log2(2*length(R))));
num_zeros = (fft_points - length(R));
Z = [Z;zeros(num_zeros,1)];
R = [R;zeros(num_zeros,1)];

n = fft_points;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ Noise estimate ----------------------------%

Z_sum=zeros(n,1);                                %  
ZZpre=zeros(n,1);                              %
Z_pre=zeros(n,1);                               %
tt=Ns:Ne;
[en,v]=dpss(Ne-Ns+1,NW,nTap);    % generate Slepian Tapers
% figure(20);
% plot(en);
en=en*sqrt(Ne-Ns+1);
ZZpre(tt)=Z(tt);                              %
for kk=1:nTap  %-------------------------%       %
    ZZpre(tt)=ZZpre(tt).*en(:,kk);       %       %
    Z_pre=fft(ZZpre)/sqrt(n);           %       %
    Z_pre2=conj(Z_pre).*Z_pre;           %       %
    Z_sum=Z_sum+real(Z_pre2);      %       %
end  %-----------------------------------%       %
SZ=Z_sum;                                  % a spectrum estimate of the pre-event noise on the vertical component

% figure(10);clf
% plot(real(ifft(SZ)),'k');
% title('pre-signal noise')
% set(gca,'XLim',[0 original_length])
% pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fourier Transform -> frequency domain
Zw = fft(Z);
Rw = fft(R);

% Do deconvolution
numerator = Rw.*conj(Zw);
denominator = Zw.*conj(Zw)+SZ;
Receiver_Fnctn_w = numerator./denominator;

% figure(20);clf
% plot(real(ifft(numerator)),'k'); hold on
% plot(real(ifft(Zw.*conj(Zw))),'b'); 
% plot(real(ifft(denominator)), 'r');
% set(gca,'XLim',[0 original_length])
% pause; 

% figure(30);clf
% plot(real(ifft(Receiver_Fnctn_w)), 'k');hold on

% multiply by delay operator to get time_shift (i.e., linear phase shift). 
n_pnts = length(Receiver_Fnctn_w);
clear i
GF = exp((i*2*pi*(time_shift*sampling_rate)/n_pnts)*[0:n_pnts-1])';
Receiver_Fnctn_w = Receiver_Fnctn_w .*GF ;

% plot(real(ifft(Receiver_Fnctn_w)), 'b');
% set(gca,'XLim',[0 original_length])
% pause; 

% Inverse Fourier Transform -> Back to the time domain
Receiver_Fnctn_t = real(ifft(Receiver_Fnctn_w));
rftn = Receiver_Fnctn_t(1:original_length);

%%%Done.































