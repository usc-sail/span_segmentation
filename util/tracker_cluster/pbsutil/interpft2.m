function im = interpft2(im,N1,N2);
% function interpft2(im,N1,N2)
%
% Function does sinc-interpolated upsampling on 2D image data using
% FFT/IFFT with zero-padding.
% The output has N1 samples in direction 1, and
% N2 samples in direction 2.
%
% Erik Bresch
% USC SPAN Group 2006.

im=interpft(im,N1,1);
im=interpft(im,N2,2);
return