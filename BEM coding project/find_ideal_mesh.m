function Ns = find_ideal_mesh(V, h_ideale)
% da V = [v1, v2, v3...] e h_ideale trovo il numero di  nodi per lato
    [lati, ~] = lati_dir(V);
    Ns = zeros(1,length(V));
    for i = 1:length(V)
        Ns(i) = ceil(lati(i)/h_ideale);
    end
end