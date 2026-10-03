clear
close all
load("finite_differences_results.mat")
load("spectral_ritz_method_results_s6.mat")



% per spectral ritz method ci servono
% eigenvalues
% eigenfunctions


for l = 1:(N+1)^2
    fig = figure();
    fig.Theme = "light";
    fcontour(eigenfunctions(l), [-1 1 -1 1])
    xlabel('x')
    ylabel('y')
    title(sprintf('Autofunzione relativa ad autovalore %f', eigenvalues(l)))
    colorbar
    axis square 
    exportgraphics(gca, "ritz6/lambda" + replace(string(eigenvalues(l)),'.',',') + ".png", 'Resolution', 300);
end

% riordino i valori in base alle autofunzioni dominanti

[valori_massimi, indici_colonna] = max(eigenvectors, [], 1);

[~, nuovo_ordine] = sort(indici_colonna);
eigenvector_ordinati = eigenvectors(:, nuovo_ordine);

eigenvalues_ordinati = eigenvalues(nuovo_ordine,:);
eigenfunctions_ordinati = eigenfunctions(nuovo_ordine,:);


fig = figure();
fig.Theme = "light";


for l = 1:(N+1)^2
    subplot(N+1,N+1,l)
    fcontour(eigenfunctions_ordinati(l), [-1 1 -1 1], 'k-')
    grid off 
    % title(sprintf('Autofunzione relativa ad autovalore %f', eigenvalues(l)))
    % colorbar
    axis square 

end



% per finite differences ci servono
% V
% lambda
% eventualmente la roba per riordinarli quindi p

% for i=4:M
%     fig = figure();
%     fig.Theme = "light";
%     contour(x_vals,x_vals,reshape(V(:,p(i)), n-2,n-2), [0,0])
%     xlabel('x')
%     ylabel('y')
%     axis equal
%     title(sprintf('Autofunzione relativa ad autovalore %f', lambda_sorted(i)))
%     colorbar
%     exportgraphics(gca, "finite_differences/lambda" + replace(string(lambda_sorted(i)),'.',',') + ".png", 'Resolution', 300);
% 
% end



