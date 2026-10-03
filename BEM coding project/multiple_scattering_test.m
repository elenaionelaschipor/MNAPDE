clear
close all

% Caso con più poligoni

% V_collection = {[0,1,1i], [1, 1.5-0.5i, 2],[2,3,3+1i], [3+1i, 3.5+1.5i, 3+2i],[3+2i,3+3i,2+3i],[2+3i, 1.5+3.5i, 1+3i],[1+3i,3i,2i], [2i, -0.5+1.5i, 1i]};

V_collection = {[0,1,1i],[2,3,3+1i]};

a = -1;
b = 4;
n_points = 150;

% V_collection = {[0,1,1i], [2,3,3+1i], [3+2i,3+3i,2+3i],[1+3i,3i,2i]};5





k = 20;
h_ideale = pi/(5*k); %così evito di non beccare le oscillazioni che sono legate a k (remark 5.9) ;

q_quadratura = ceil(3*k); % idem

[Nodes, Midpoints, Hj, Dirs] = create_nodes(V_collection, h_ideale);

d = [1/2, -sqrt(3)/2];


plane_wave = @(z) cos(k*(real(z)*d(1)+imag(z)*d(2))) + 1i*sin(k*(real(z)*d(1)+imag(z)*d(2))); % exp(1i*k*(real(z)*d(1) + imag(z)*d(2)));
x0 = 1.5+1.5i;
sol_fondamentale_x0 = @(z) 0.25i*besselh(0,k*abs(z-x0));


onda_incidente = plane_wave;

u =build_2D_solution(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);

x = linspace(a,b,n_points);
[xx,yy] = meshgrid(x,x);               
zz = xx + 1i * yy; 
u_inc =  arrayfun(@(z) onda_incidente(z), zz) ;
u_tot = u + u_inc;



fig = figure();
fig.Theme = "light";
subplot(3,3,1)
imagesc(x,x,real(u_inc))
title("R(u_{inc})")
set(gca, 'YDir', 'normal')   
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
subplot(3,3,2)
imagesc(x,x,imag(u_inc))
title("I(u_{inc})")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

subplot(3,3,3)
imagesc(x,x,abs(u_inc))
title("|u_{inc}|")
set(gca, 'YDir', 'normal') 
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)


subplot(3,3,4)
imagesc(x,x,real(u))
title("R(u_{scat})")
set(gca, 'YDir', 'normal') 
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

subplot(3,3,5)
imagesc(x,x,imag(u))
title("I(u_{scat})")
set(gca, 'YDir', 'normal')   
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

subplot(3,3,6)
imagesc(x,x,abs(u))
title("|u_{scat}|")
set(gca, 'YDir', 'normal') 
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

subplot(3,3,7)
imagesc(x,x,real(u_tot))
title("R(u_{tot})")
set(gca, 'YDir', 'normal')   
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

subplot(3,3,8)
imagesc(x,x,imag(u_tot))
title("I(u_{tot})")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

subplot(3,3,9)
imagesc(x,x,abs(u_tot))
title("|u_{tot}|")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)

exportgraphics(fig, "multiple_scattering.png", 'Resolution', 300);

% per le animazioni
% MyGiffer(xx,yy,u,V_collection,"u_scatt_sol_fond.gif")
% MyGiffer(xx,yy,u_inc,V_collection, "u_inc_sol_fond.gif")
% MyGiffer(xx,yy,u_tot, V_collection,"u_tot_sol_fond.gif")

