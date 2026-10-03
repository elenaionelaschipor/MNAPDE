function [A, F] = build_collocation(V, h_ideale, k, dato_bordo, q_quadratura)
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

            A(j,m) = 0.25i*hm*sum(w.*besselh(0, k*abs(pm + (   hm*nodi_quadratura   )*tau_m - xj)));

        else
            hj = Hj(j);
            eugamma = -psi(1);
            A(j,j) = (0.5i - (1/pi)*eugamma + (1/pi)*log(2))*hj/2 - (1/pi)*(hj/2)*log(k*hj/2) + (1/pi)*hj/2;
            % A(j,j) = 0.25i*hj - 1/(2*pi)*hj*log(k*hj/2) + + (1/pi)*hj/2;
            % A(j,j) = - 1/(2*pi)*hj*log(k*hj/2) + (1/pi)*hj/2;
        end 
    end
    F(j) = dato_bordo(xj);
end


end