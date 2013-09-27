% comparison of TV projection algorithms
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

name = 'lena';
n = 256;

rep = 'results/tv-proj-comp/';
if not(exist(rep))
   mkdir(rep);
end

M = load_image(name, n);
M = rescale(M);

tv = compute_total_variation(M);
tvtgt = tv/2;

niter = 400;
niter = 2000;


%% pre-compute a solution with large # iterations
clear options;
options.niter = niter;
options.tol = 0;
options.nesterov = 1;
[Minf,err,tverr0,u] = perform_tv_projection_fb(M, tvtgt, options); 
tverr0 = tverr0 - tvtgt;

niter1 = niter/2;
%% Forward Backward %%
options.Mtgt = Minf;
options.niter = niter1;
options.nesterov = 0;
[M0,err,tverr0,u, errtgt0] = perform_tv_projection_fb(M, tvtgt, options); 
tverr0 = tverr0 - tvtgt;
errtgt0 = errtgt0 / norm(M, 'fro');

%% Nesterov %%
options.niter = niter1/2;
options.nesterov = 1;
[M2,err,tverr2,u, errtgt2] = perform_tv_projection_fb(M, tvtgt, options); 
tverr2 = tverr2 - tvtgt;
errtgt2 = errtgt2 / norm(M, 'fro');

%% Combettes %%
options.niter = niter1;
options.xtgt = Minf;
[M1,tverr1, err, errtgt1] = perform_tv_projection(M, tvtgt, options);
errtgt1 = errtgt1 / norm(M, 'fro');

%% display decay of TV norm
nesterovx = 2:2:niter1;
clf;
hold on;
h = plot( log10(tverr0/tvtgt), 'k-' );
set(h, 'LineWidth', 3);
h = plot( log10(tverr1/tvtgt), 'k--' );
set(h, 'LineWidth', 3);
h = plot( nesterovx, log10(tverr2/tvtgt), 'k-.' );
set(h, 'LineWidth', 3);
axis( [1 niter1 min(log10( tverr2/tvtgt ))  max(log10( tverr2/tvtgt )) ] );
legend('Forward-Backward', 'Subgradient proj.', 'Nesterov');
saveas(gcf, [rep name '-tv-error.eps'], 'eps');

%% display decay of error
clf;
hold on;
h = plot( log10( errtgt0 ), 'k-' );
set(h, 'LineWidth', 3);
h = plot( log10( errtgt1 ), 'k--' );
set(h, 'LineWidth', 3);
h = plot( nesterovx, log10( errtgt2), 'k-.' );
set(h, 'LineWidth', 3);
axis( [1 niter1 min(log10( errtgt2 ))  max(log10( errtgt2 )) ] );
legend('Forward-Backward', 'Subgradient proj.', 'Nesterov');
saveas(gcf, [rep name '-l2-error.eps'], 'eps');