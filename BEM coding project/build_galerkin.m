function [A, F] = build_galerkin(V, h_ideale, k, dato_bordo, q_quadratura)
%% V = vertici del poligono 
% h_ideale = taglia (circa) della mesh
% k = k della eq di helmholtz
% dato_bordo = dato al bordo

% Ns = find_ideal_mesh(V, h_ideale); 
[Nodes, Midpoints, Hj, Dirs] = create_nodes(V, h_ideale);  % ok

N = length(Nodes);
% disp("N è" + string(N))
% disp("M è" + string(M))
A = zeros(N);
[nodi_quadratura,w] = gaussquad(q_quadratura);
F = zeros(N,1);
for j=1:N

    xj = Midpoints(j);
    for m=1:N
        if j ~=m
            hm = Hj(m);
            pm = Nodes(m);
            tau_m = Dirs(m);
            hj = Hj(j);
            pj = Nodes(j);
            tau_j = Dirs(j);
            

            A(j,m) = 0.25i*hm*sum(w'.*( hj*sum( w.*besselh(0, k*abs(pm + (   hm*nodi_quadratura   )*tau_m -pj - (   hj*nodi_quadratura'   )*tau_j )))));

        else
            hj = Hj(j);
            tau_j = Dirs(j);
            % 
            % eugamma = -psi(1);
            % A(j,j) = 0.25i*(1+ (2/pi)*eugamma)*hj^2 - hj^2/(2*pi)*log(k*abs(tau_j)/2) ...
            %    - hj^2*(log(hj)-3/2);

            A(j,j) = 1i*hj^2/2 * sum(w.*(1-nodi_quadratura).*besselh(0, k*hj*nodi_quadratura) );

        end 
    end
    F(j) = hj*sum(w.*dato_bordo ( (pj + (   hj*nodi_quadratura   )*tau_j   )));
end


end