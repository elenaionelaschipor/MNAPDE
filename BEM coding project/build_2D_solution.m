function u = build_2D_solution(V, h_ideale, k, onda_incidente,a,b,n_points, q_quadratura)
%% onda incidente: funzione del tipo @(x) ...

% Ns = find_ideal_mesh(V, h_ideale);
[Nodes, Midpoints, Hj, Dirs] = create_nodes(V, h_ideale);

gD = @(x) - onda_incidente(x);        % - traccia di onda incidente 

[A, F] = build_collocation(V, h_ideale, k, gD, q_quadratura);

PSI = A\F;  % soluzione sul bordo


x = linspace(a,b,n_points);
[xx,yy] = meshgrid(x,x);                % Matrices of x and y coordinates where to evaluate and plot field
zz = xx + 1i * yy; 
[s,w] = gaussquad(q_quadratura);

Z = reshape(zz, [], 1);
Y = Nodes-Dirs.*Hj.*s;
Y = reshape(Y, 1, size(Y,1), size(Y,2));

dist = abs(Z-Y);

G = besselh(0,1, k*dist);
SLP_VALS= G.*transpose(w);
cici =sum(SLP_VALS, 2);
cici = reshape(cici, n_points^2, []);
P = Hj .* cici;



u = 0.25i * sum(transpose(PSI).*P, 2);   % somma su M

u = reshape(u, n_points, n_points);

% P = Hj.*sum(w.*besselh(0,1,k*abs(Z - Nodes-Dirs.*Hj.*s)));



%u = 0.25i*arrayfun(@(z) sum(transpose(PSI).*(SLP(z, Nodes, Dirs, Hj, s,w, k))), zz);

end


