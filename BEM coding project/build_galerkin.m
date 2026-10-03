function [A, F] = build_galerkin(V, h_ideale, k, dato_bordo, q_quadratura)
%% V = vertici del poligono 
% h_ideale = taglia (circa) della mesh
% k = k della eq di helmholtz
% dato_bordo = dato al bordo


[Nodes, Midpoints, Hj, Dirs] = create_nodes(V, h_ideale); 

N = length(Nodes);
A = zeros(N);
[nodi_quadratura,w] = gaussquad(q_quadratura);
F = zeros(N,1);
for j=1:N
    for m=1:N
        if j ~=m  % caso facile, faccio integrali e basta
            hm = Hj(m);
            pm = Nodes(m);
            tau_m = Dirs(m);
            hj = Hj(j);
            pj = Nodes(j);
            tau_j = Dirs(j);
            
            A(j,m) = 0.25i*hm*sum(w'.*( hj*sum( w.*besselh(0, k*abs(pm + (   hm*nodi_quadratura   )*tau_m -pj - (   hj*nodi_quadratura'   )*tau_j )))));

        else % uso formula in 5.2.1 delle dispense
            hj = Hj(j);
            tau_j = Dirs(j);
            A(j,j) = 1i*hj^2/2 * sum(w.*(1-nodi_quadratura).*besselh(0, k*hj*nodi_quadratura) );
        end 
    end
    F(j) = hj*sum(w.*dato_bordo ( (pj + (   hj*nodi_quadratura   )*tau_j   )));
end
end