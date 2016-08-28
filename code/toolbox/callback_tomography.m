function y = callback_tomography(x, dir, options)

% callback_tomography - callback for tomography
%
%   options.mask gives the fourier sub-sampling mask.
%
%   y = callback_tomography(x, dir, options);
%
%   Copyright (c) Gabriel Peyre

mask = getoptions(options, 'mask', 1, 1);
n = size(x,1);

if dir>0
    y = fft2(x) / n; y(mask==0) = 0;
else
    x(mask==0) = 0; 
    y = real(ifft2(x)) * n;
end