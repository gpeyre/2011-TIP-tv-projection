% test for the resolution of inverse problems using projection
%
%	Copyright (c) 2009 Gabriel Peyre
%
%	For more information, see
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyre,
%		Proceedings of ICIP'09, Nov. 2009.
%	and
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyré,
%		Preprint Hal-00380491, May 2009.


path(path, 'toolbox/');

problem = 'deblurring';
problem = 'inpainting';
problem = 'tomography';

name = 'boat';
name = 'cameraman';
name = 'lena';
name = 'phantom';
n = 256;

rep = 'results/tv-proj-ip/';
if not(exist(rep))
    mkdir(rep);
end

n0 = [];
if strcmp(name, 'phantom')
    n0 = n;
end

%% load image
M0 = load_image(name,n0);
M0 = crop(M0, n);
M0 = rescale(M0);

if strcmp(name, 'cameraman')
    % remove first row
    s= 2;
    M0(1:s,:) = M0(2*s:-1:s+1,:);
end

%% set up IP
switch problem
    case 'inpainting'
        % create random mask
        rho = .7; % ratio of missing pixels
        mask = zeros(n);
        sel = randperm(n*n); sel = sel(1:round(rho*end));
        mask(sel) = 1;
        options.mask = mask;
        mu = 1;
        cbk = @callback_inpainting;
        sigma = .05; % no noise
        taumult = .6;
    case 'deblurring'
        % create blurring mask
        tau = 4;
        x = [0:n/2-1 -n/2:-1];
        [Y,X] = meshgrid(x,x);
        filter = exp( -(X.^2+Y.^2)/(2*tau^2) );
        filter = filter / sum(filter(:));
        options.filter = fft2(filter);
        mu = 1.9;
        cbk = @callback_convolution;
        sigma = .02; % noise
        taumult = .3;
    case 'tomography'
        % create the sub-sampling mask
        nrays = 22;
        Theta = linspace(0,pi,nrays+1); Theta(end) = [];
        mask = zeros(n);
        for theta = Theta
            t = linspace(-1,1,3*n)*n;
            x = round(t.*cos(theta)) + n/2+1; y = round(t.*sin(theta)) + n/2+1;
            I = find(x>0 & x<=n & y>0 & y<=n); x = x(I); y = y(I);
            mask(x+(y-1)*n) = 1;
        end
        options.mask = fftshift(mask);
        %
        mu = 1.9;
        cbk = @callback_tomography;
        sigma = .02; % noise
        taumult = .3; % target TV
end

%% set ip observations
y = feval(cbk, M0, +1, options);
y = y + randn(size(y))*sigma;

tv0 = compute_total_variation(M0);
tv = tv0 * taumult;

% iterations for ip
niter = 200;
% iterations for projection
options.niter = 30;
options.verb = 0;
options.tol = 1e-6;
% initialization
M1 = zeros(n);
options.u = [];
err2 = [];
Msvg = [];
for i=1:niter
    progressbar(i,niter);
    % gradient step   M1 <-  M1 + mu*Astar( y-A*M1 )
    r = y - feval(cbk, M1, +1, options);
    err2(end+1) = norm(r, 'fro');
    M1 = M1 + mu * feval(cbk, r, -1, options);
    % projection step
    [M1,err,tverr,u] = perform_tv_projection_fb(M1, tv, options); 
    M1 = clamp(M1);
    clf; imageplot(M1); drawnow;
    Msvg(:,:,i) = M1;
end

options.base_str = [rep name '-' problem '-'];
save_image(clamp({M0 y M1}), {'original' 'observation' 'solution'}, options);

clf;
h = plot(err2, 'k');
set(h, 'LineWidth', 2);
axis tight;
saveas(gcf, [rep name '-' problem '-energy.eps'], 'eps');

% plot decay of the energy
m = niter/2;
err = sum(sum( (Msvg(:,:,1:m)-repmat(Msvg(:,:,end), [1 1 m])).^2 ));
err = sqrt(err(:)) / norm(M0, 'fro');
h = plot(log10(err), 'k');
set(h, 'LineWidth', 2);
axis tight;
saveas(gcf, [rep name '-' problem '-error.eps'], 'eps');
