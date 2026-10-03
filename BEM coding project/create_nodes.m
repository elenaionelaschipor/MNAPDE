function [Nodes, Midpoints, Hj, Dirs] = create_nodes(V_collection, h_ideale)
%% V_collection = cell con array che sono i vertici dei poligoni
% h ideale = numero 
    Ns_totali = cell(1,length(V_collection));
    N = 0;
    ind = 1;

    for V_ = V_collection  % per ogni poligono...
        V = V_{:};
        Ns = find_ideal_mesh(V, h_ideale);  % trovo mesh con h_ideale, Ns = [a,b,c,...] contiene il numero di nodi per ogni lato
        Ns_totali(ind) = {Ns};
        ind = ind + 1; 
        N = N + sum(Ns);  % N sarà il numero totale di nodi
    end
       
    
    Nodes = zeros(1,N);  % pj
    Midpoints = Nodes; % xj
    Dirs = Nodes; %tau_j
    Hj = Nodes; %h_j

    ind = 1;
    idx = 1;

    for V_ = V_collection  % per ogni poligono...
        V = V_{:};
        Ns = Ns_totali{ind};
        [lati, direzioni] = lati_dir(V);
        
         
        Nodes(idx) = V(1); % per ogni lato...
        for j = 1:length(lati)
            N = Ns(j);
            h = lati(j)/N;
            tau = direzioni(j);
            for i = 0:(N-1)
                Nodes(idx) = V(j) + h*tau*i;
                Midpoints(idx) = Nodes(idx) + h/2*tau;
                Hj(idx) = h;
                Dirs(idx) = tau;        
    
                idx = idx+1;
            end  
            
        end
        ind = ind+1;
    end
end