clear
close all
load("spectral_ritz_method_results_s6.mat")  % workspace ottenuto dopo spectral_ritz_method.m



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



% load("finite_differences_results.mat")  % workspace ottenuto dopo finite_differences.m
 
% per finite differences ci servono
% V
% lambda
% eventualmente per riordinarli p

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



