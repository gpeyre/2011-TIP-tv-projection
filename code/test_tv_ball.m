%%
% Display the 4D TV Ball in 3D.

path(path, 'toolbox/');

% dimension 2/3
d = 3;
% type 1=TV, 2=Sobolev, 3=tviso
t = 1;
t = 3;

randn('state', 12341);
U = randn(d+1);
U(:,1) = 1;
[U,R] = qr(U);
U = U(:,2:end);

if d==2
    %%%%%% 3D %%%%%%%%
    
    n = 512;
    x = linspace(-1,1,n);
    [A,B] = ndgrid(x,x);
    V = zeros(n,n,3);
    for i=1:3
        V(:,:,i) = U(i,1)*A + U(i,2)*B;
    end
    
    T = zeros(n,n);
    for i=1:3
        j = mod(i,3)+1;
        T = T + abs(V(:,:,i) - V(:,:,j)).^t;
    end
    
    
    clf;
    imageplot(T<=2.2);
    
else
    %%%%%% 3D %%%%%%%%
    
    
    
    n = 105;
    x = linspace(-1,1,n);
    [A,B,C] = ndgrid(x,x,x);
    V = zeros(n,n,n,4);
    for i=1:4
        V(:,:,:,i) = U(i,1)*A + U(i,2)*B + U(i,3)*C;
    end
    
    T = zeros(n,n,n);
    
    if t==1 || t==2
        for i=1:4
            j = mod(i,4)+1;
            T = T + abs(V(:,:,:,i) - V(:,:,:,j)).^t;
        end
    elseif t==3
        % 1 2 
        % 3 4
        V = reshape(V, n,n,n,2,2);
        for i=1:2
            for j=1:2
                i1 = mod(i,2)+1;
                j1 = mod(j,2)+1;
                T = T + sqrt( (V(:,:,:,i,j) - V(:,:,:,i,j1)).^2 + (V(:,:,:,i,j) - V(:,:,:,i1,j)).^2 );
            end
        end        
    end
    
    if t==1 || t==3
        tau = 2.1;
    else
        tau = 2;
    end
    clf;
    plot_isosurface(T,tau);
    view(-72,16);
    
end
saveas(gcf, ['tv-ball-' num2str(d) 'd-t' num2str(t) '.png'], 'png');