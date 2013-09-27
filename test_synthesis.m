% test for texture synthesis using projected noise
%
%	Copyright (c) 2009 Gabriel Peyre
%
%	For more information, see
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyre
%		Proceedings of ICIP'09, Nov. 2009.
%	and
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyre,
%		Preprint Hal-00380491, May 2009.

histo = 'linear';
n = 256;

path(path, 'toolbox/');

rep = 'results/tv-synthesis/';
if not(exist(rep))
    mkdir(rep);
end

% type of spacial constraint
constraint = 'hist';
constraint = 'std';

% projection on unit ball
proj = @(x)(x-mean(x(:)))/norm(x-mean(x(:)),'fro');

clear options;
niter_synth = 14;

randn('state', 1234556);
M0 = randn(n);
if strcmp(constraint, 'hist')
    M0 = perform_histogram_equalization( M0, histo);
else
    M0 = proj(M0);
end
tv0 = compute_total_variation(M0);

tvlist = tv0 ./ [4 8 16 32 64];
tvlist = tv0 ./ [4 8 16 32];
niterlist = [100 200 400 800 1600]*4;
options.base_str = [rep 'tv-synthesis-' constraint '-'];

a = 1.2/n;
Msvg = {};
for k=1:length(tvlist)
    options.niter = niterlist(k);
    tvtgt = tvlist(k);    
    M1 = M0;
    for i=1:niter_synth
        tvtgt = tvlist(k);
        [M1,err,tverr,options.u] = perform_tv_projection_fb(M1, tvtgt, options);
        tverr = tverr - tvtgt;
        % equalize
        if strcmp(constraint, 'hist')
            M1 = perform_histogram_equalization( M1, histo );
        else
            M1 = proj(M1);
        end
    end
    %
    Msvg{end+1} = M1;
    save_image( rescale(clamp(M1,-a,a)), num2str(round(tv0 / tvtgt)), options);
end