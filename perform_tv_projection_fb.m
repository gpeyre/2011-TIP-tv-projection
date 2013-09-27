function [M1,err,tv, u, errtgt] = perform_tv_projection_fb(M, tv0, options)

% perform_tv_projection_fb - compute projection on the TV ball
%
%   [M1,err,tv, u] = perform_tv_projection_fb(M, tv0, options);
%
%   Solve 
%       min_M1  |M-M1|  subj.to  TV(M1)<tv0
%
%   using an iterative soft threshoding on the dual.
%
%   The number of iterations is options.niter
%   You can give an initial vector field in options.u.
%
%   You can set options.nesterov=1 to use a faster Nesterov algorithm.
%
%   Copyright (c) 2009 Gabriel Peyre
%
%	For more information, see
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyré,
%		Proceedings of ICIP'09, Nov. 2009.
%	and
%		Total Variation Projection with First Order Schemes
%		Jalal Fadili and Gabriel Peyré,
%		Preprint Hal-00380491, May 2009.

options.null = 0;

niter = getoptions(options, 'niter', 1000);
gradop = getoptions(options, 'gradop', @grad);
divop = getoptions(options, 'divop', @div);
mu = getoptions(options, 'mu', 0.25);
disp = getoptions(options, 'display', 1);
ndisp = getoptions(options, 'ndisp', 100);
u = getoptions(options, 'u', []);
Mtgt = getoptions(options, 'Mtgt', []);
verb = getoptions(options, 'verb', 1);
tol = getoptions(options, 'tol', 1e-6);
nesterov = getoptions(options, 'nesterov', 1);


n = size(M,1);

if isempty(u)
    u = zeros(n,n,2);
end
if nesterov
    u0 = zeros(n,n,2);
    g = zeros(n,n,2);
    A = 0;
end

err = [];
tv = [];
nM = norm(M,'fro');

for i=1:niter
    if verb
        progressbar(i,niter);
    end
    
    if nesterov
        
        a = ( mu + sqrt(mu^2 + 4*mu*A) )/2;
        
        % prox_{A}(u0 - g)
        v = u0-g;
        v = v - perform_l1ball_projection( v, tv0*A );
        %
        y = (A*u + a*v)/(A+a);
        % prox_{mu/2}( y-mu/2*grad(M-div(y)) )
        u = y - mu/2 * grad(M-div(y,options),options);
        u = u - perform_l1ball_projection( u, tv0*mu/2 );
        %
        Mnew = M-div(u,options);
        g = g + a * grad(Mnew);
        if i>1 && norm(M1-Mnew, 'fro')/nM<tol
            M1 = Mnew; 
            progressbar(niter,niter);
            break;
        end
        M1 = Mnew;
        
        A = A + a;
    
    else
        
        % descent step
        Mnew =  M - feval( divop, u, options);
        if i>1 && norm(M1-Mnew, 'fro')/nM<tol
            M1 = Mnew; break;
        end
        M1 = Mnew;
        u = u - mu*feval( gradop, M1, options );
        % prox step
        u = u - perform_l1ball_projection( u, tv0*mu );
       
        
    end
    
    
    % display
    if disp && mod(i, max(round(niter/ndisp),1) )==1
        imageplot(M1); drawnow;
    end
    % monitor error
    err(i) = 1/2*norm(M1,'fro')^2;  % + tv0 * max(max( sqrt(sum(u.^2,3)) ));
    tv(i) = compute_total_variation(M1);
    if not(isempty(Mtgt))
        errtgt(i) = norm(M1-Mtgt, 'fro');
    end
end