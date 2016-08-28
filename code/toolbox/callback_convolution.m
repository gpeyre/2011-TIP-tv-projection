function y = callback_convolution(x, dir, options)

% callback_convolution - callback for convolution
%
%   options.filter gives the FOURIER transform of the filter, 
% that should be symmetric.
%
%   Copyright (c) Gabriel Peyre

filter = getoptions(options, 'filter', 1, 1);
y = real( ifft2(fft2(x).*filter) );