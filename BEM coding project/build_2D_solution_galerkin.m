function u = build_2D_solution_galerkin(V, h_ideale, k, onda_incidente,a,b,n_points, q_quadratura)
%% onda incidente: funzione del tipo @(x) ...

% Ns = find_ideal_mesh(V, h_ideale);
[Nodes, Midpoints, Hj, Dirs] = create_nodes(V, h_ideale);

gD = @(x) - onda_incidente(x);  % - traccia di onda incidente 

[A, F] = build_galerkin(V, h_ideale, k, gD, q_quadratura);

PSI = A\F;  % soluzione sul bordo

x = linspace(a,b,n_points);
[xx,yy] = meshgrid(x,x);               
zz = xx + 1i * yy; 
[s,w] = gaussquad(q_quadratura);

Z = reshape(zz, [], 1);
Y = Nodes-Dirs.*Hj.*s;
Y = reshape(Y, 1, size(Y,1), size(Y,2));

% formula di rappresentazione
dist = abs(Z-Y); 
G = besselh(0,1, k*dist);
SLP_VALS= G.*transpose(w);
c =sum(SLP_VALS, 2);
c = reshape(c, n_points^2, []);
P = Hj .* c;

u = 0.25i * sum(transpose(PSI).*P, 2);   % somma su M

u = reshape(u, n_points, n_points);

end


