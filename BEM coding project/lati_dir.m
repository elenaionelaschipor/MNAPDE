function [lati, direzioni] = lati_dir(V)
lati = zeros(1, length(V));
direzioni = lati;
lati(end) = abs(V(end) - V(1));
direzioni(end) = (V(1) - V(end))/abs(V(end) - V(1));
for i = 2:length(V)
    lati(i-1) = abs(V(i) - V(i-1));
    direzioni(i-1) = ( V(i) - V(i-1) )/abs(V(i) - V(i-1));
end
end
