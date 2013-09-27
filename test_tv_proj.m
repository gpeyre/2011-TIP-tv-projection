% Test for total variation projection.
%
%	Copyright (c) 2009 Gabriel Peyre
%
%	For more information, see
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyré,
%		Proceedings of ICIP'09, Nov. 2009.
%	and
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyré,
%		Preprint Hal-00380491, May 2009.

path(path, 'toolbox/');

name = 'lena';
n = 256;

rep = 'results/tv-proj/';
if not(exist(rep))
    mkdir(rep);
end

M = load_image(name, n);
M0 = rescale(M,.05,.95);
sigma = .0;
M = M0 + randn(n)*sigma;

options.base_str = [rep name '-'];
save_image(clamp({M M0}), {'noisy' 'original'}, options);

tv0 = compute_total_variation(M0);
tv = compute_total_variation(M);

tv_list = tv0 ./ [.5 2 4 8 16];
niter = [200 400 1000 4000]*2;
u = [];
tverr = {}; err = {};
clf;
for k=1:length(tv_list);
    tvtgt = tv_list(k);
    options.niter = niter(k);
    options.u = u;
    [M1,err{k},tverr{k},u] = perform_tv_projection_fb(M, tvtgt, options); 
    tverr{k} = tverr{k} - tvtgt;
    imageplot(clamp(M1), num2str(tv0 / tvtgt), 2, 2, min(k,4) );
    save_image(clamp(M1), num2str(round(tv0 / tvtgt)), options);
end