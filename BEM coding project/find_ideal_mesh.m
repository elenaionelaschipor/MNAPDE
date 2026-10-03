function Ns = find_ideal_mesh(V, h_ideale)
    [lati, ~] = lati_dir(V);
    Ns = zeros(1,length(V));
    for i = 1:length(V)
        Ns(i) = ceil(lati(i)/h_ideale);
    end
end