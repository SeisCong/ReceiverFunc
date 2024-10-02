function [X] = gaussf(Y,sampling_rate,a)

% function [X] = gaussf(Y,sampling_rate,a)
%
% This is a gaussian low pass filter.
% 
%             Y: data to be filtered.
%                data down columns.
% sampling_rate: sampling rate in samples/s
%             a: gaussian parameter.
%                controls corner.
%
% A good rule of thumb for a is that the
% approximate corner is half the value for a:
%
% 'a' value    |  approximate corner
% ________________________________
%    10               4.8
%    5                2.4
%    2.5              1.2
%    1.25             .6
%    .625             .3
%    .5               .24
%    .4               .2
%    .2               .1

%copyright Jason Crosswhite (c) 2002
%jason@newberry.uoregon.edu

%if row vector -> column vector
%shouldn't have to worry about things
%like this, as seismograms should always
%be down columns, but I kept forgetting
%with a vector of size 1xn or nx1, so
%I put this in here. 
if (size(Y,1) == 1)
     Y = Y(:);
end
original_points = size(Y,1);

%get nearest ^2 for fft.
points = 2^(ceil(log2(original_points)));  

%into the frequency domain
w_Y = fft(Y,points);
nyquist = sampling_rate/2;
dw = (2*nyquist)/(points-1);
freqs = [-nyquist:dw:nyquist]*2*pi;
gaussian = exp(-(freqs.^2)/(4*a^2));

%need to normalize the gaussian to unit amplitude.
%constant = max(real(ifft(gaussian)));

%make a matrix full of gaussians down columns
mat = repmat(fftshift((gaussian)'),1,size(Y,2)); 

%filter data.
w_X = w_Y.*mat;

%back to time domain
X = real(ifft(w_X));
X = X(1:original_points,:);
