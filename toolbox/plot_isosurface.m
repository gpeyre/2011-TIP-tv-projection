function plot_isosurface(M,rho)

if nargin<2
    rho = .5;
end

F = isosurface( M,rho );
p = patch(F);
isonormals( M,p );
set(p, 'FaceColor', 'red', 'EdgeColor', 'none'); 

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'ZTick', []);
lighting phong;
% alpha(0.7);    
camproj('perspective');
view(3);
% axis equal;
% axis([1 size(W,1) 1 size(W,2) 1 size(W,3)]);
% daspect([1 1 1]);
% cameramenu;
axis equal; axis off;
camlight;