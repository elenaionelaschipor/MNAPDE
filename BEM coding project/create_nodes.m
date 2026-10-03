function [Nodes, Midpoints, Hj, Dirs] = create_nodes(V_collection, h_ideale)
    % figure()
    % plot_polygon(V)
    % hold on

    Ns_totali = cell(1,length(V_collection));
    N = 0;
    ind = 1;
    for V_ = V_collection
        V = V_{:};
        Ns = find_ideal_mesh(V, h_ideale);
        Ns_totali(ind) = {Ns};
        ind = ind + 1; 
        N = N + sum(Ns);
    end
        

    
    Nodes = zeros(1,N);
    Midpoints = Nodes;
    Dirs = Nodes;
    Hj = Nodes;

    ind = 1;
    idx = 1;
    for V_ = V_collection
        V = V_{:};
        Ns = Ns_totali{ind};
        [lati, direzioni] = lati_dir(V);
        
        
        Nodes(idx) = V(1);
        for j = 1:length(lati)
            % disp("lato" + string(j))
            N = Ns(j);
            h = lati(j)/N;
            tau = direzioni(j);
            for i = 0:(N-1)
                % disp(idx)
                Nodes(idx) = V(j) + h*tau*i;
                Midpoints(idx) = Nodes(idx) + h/2*tau;
                Hj(idx) = h;
                Dirs(idx) = tau;        
    
                idx = idx+1;
                 
                % plot(real(Nodes(1:idx-2)), imag(Nodes(1:idx-2)), 'or', MarkerFaceColor='r')
                % plot(real(Midpoints(1:idx-2)), imag(Midpoints(1:idx-2)), 'ob', MarkerFaceColor='b')
            end  
            
        end
        ind = ind+1;
    end
    
    % V = [1, 2, 2+3i , 1i];
    % Ns = find_ideal_mesh(V, 0.1);
    % 
    % figure()
    % [Nodes, Midpoints, Hj, Dirs] = create_nodes(V, Ns);
    % plot_polygon(V)
    % % plot(V, '-or')
    % hold on
    % plot(Nodes, 'or', MarkerFaceColor='r')
    % plot(Midpoints, 'ob', MarkerFaceColor='b')
    % 

end

% 
% V_collection = {[0,1,1i],[2,2+0.3i,3+0.5i,1.5i]};
% h_ideale = 0.1;
% figure()
% [Nodes, Midpoints, Hj, Dirs] = create_nodes(V_collection, h_ideale);
% plot_polygon(V_collection)
% % plot(V, '-or')
% hold on
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, 'ob', MarkerFaceColor='b')