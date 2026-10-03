clear
close all

n = 100;
mu = 0.225;
h = 2/(n-3);
M = 90; 


G = numgrid('S', n+2);
D = delsq(G);

N = D;
L = D; 
% boundary points
bleft = G(3:n,3) ; % left
bright = G(3:n,n); % right
btop =  G(3,3:n)'; %top
bbottom = G(n,3:n)'; %bottom 

% ghost points
gleft = G(3:n,2); 
gright = G(3:n,n+1);
gtop = G(2,3:n)';
gbottom = G(n+1,3:n)';


N(bleft, bleft) = N(bleft, bleft)/2;
N(bright, bright) = N(bright, bright)/2;
N(btop, btop) = N(btop, btop)/2;
N(bbottom, bbottom) = N(bbottom, bbottom)/2;


L([gleft;gright;gtop;gbottom],:) = 0;
for j=gleft(1 :end-1)' % left boundary
    L([j,j+1] , [j,j+1,j+2*n,j+2*n+1]) = L([j,j+1] , [j,j+1,j+2*n,j+2*n+1]) + (mu-1)/2*[1,-1,-1,1;-1, 1,1,-1] ;
end

for j=gright(1 :end-1)' % right boundary
    L([j,j+1] , [j,j+1,j-2*n,j-2*n+1]) = L([j,j+1] , [j,j+1,j-2*n,j-2*n+1]) + (mu-1)/2*[1,-1,-1,1;-1, 1,1,-1] ;
end



for j=gtop(1 :end-1)' % top boundary
    L([j,j+n] , [j+n,j,j+n+2,j+2]) = L([j,j+n] , [j+n,j,j+n+2,j+2]) - (mu-1)/2*[1,-1,-1,1;-1, 1,1,-1] ;
end


for j=gbottom(1 :end-1)' % bottom boundary
    L([j,j+n] , [j+n,j,j+n-2,j-2]) = L([j,j+n] , [j+n,j,j+n-2,j-2]) - (mu-1)/2*[1,-1,-1,1;-1, 1,1,-1] ;
end



A = N*L;



A(gleft,:) = 0; % left ghost points
for i=gleft'
    A(i, [i+n,i,i+n-1,i+n+1,i+2*n] ) = [2*(1+mu), -1, -mu, -mu, -1];
end

A(gright,:) = 0; % right ghost points
for i=gright'
    A(i, [i-n,i,i-n-1,i-n+1,i-2*n] ) = [2*(1+mu), -1, -mu, -mu, -1];
end


A(gtop,:) = 0; % top ghost points
for i=gtop'
    A(i, [i+1,i,i+n+1,i-n+1,i+2] ) = [2*(1+mu), -1, -mu, -mu, -1];
end


A(gbottom,:) = 0; % bottom ghost points
for i=gbottom'
    A(i, [i-1,i,i+n-1,i-n-1,i-2] ) = [2*(1+mu), -1, -mu, -mu, -1];
end



% elimino ghost nodes

phys = G(3:n,3:n); 
phys = phys(:); % put all physical nodes in a vector
ghost = [gleft; gright; gtop; gbottom] ;
A0 = A (phys, phys) - A (phys,ghost)/A (ghost , ghost )*A(ghost , phys);



% RHS
B= h^4*speye(n^2);
B(bleft, bleft) = B(bleft,bleft)/2;
B(bright, bright) = B(bright,bright)/2;
B(btop,btop) = B(btop,btop)/2;
B(bbottom, bbottom)=B(bbottom,bbottom)/2;

B0 = B(phys, phys);

[V, lambda] = eigs(A0, B0, M, 'smallestabs');
[lambda_sorted,p] = sort(diag(lambda));

x_vals = -1:2/(n-3):1;

V = real(V);
for i=4:M
    figure()
    contour(x_vals,x_vals,reshape(V(:,p(i)), n-2,n-2), [0,0])
    axis equal
end

