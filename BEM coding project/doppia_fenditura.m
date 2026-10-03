clear
close all


a = -1.5;
b = 1.5;
n_points = 150; % ocio che se son troppi esplode

V_collection = {[-1.5,-0.30,-0.3+0.005i, -1.5+0.005i], [0,0.1,0.1+0.005i, 0.005i], [0.4,1.5,1.5+0.005i, 0.4+ 0.005i]};

k = 40;
h_ideale = pi/(5*k); %così evito di non beccare le oscillazioni che sono legate a k (remark 5.9) ;

q_quadratura = ceil(3*k); % idem


% Ns = find_ideal_mesh(V_collection, h_ideale);
[Nodes, Midpoints, Hj, Dirs] = create_nodes(V_collection, h_ideale);

d = [0.5, sqrt(3)/2];


plane_wave = @(z) cos(k*(real(z)*d(1)+imag(z)*d(2))) + 1i*sin(k*(real(z)*d(1)+imag(z)*d(2))); % exp(1i*k*(real(z)*d(1) + imag(z)*d(2)));
x0 = 1.5+1.5i;
sol_fondamentale_x0 = @(z) 0.25i*besselh(0,k*abs(z-x0));


onda_incidente = plane_wave;

u =build_2D_solution(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);

x = linspace(a,b,n_points);
[xx,yy] = meshgrid(x,x);                % Matrices of x and y coordinates where to evaluate and plot field
zz = xx + 1i * yy; 
u_inc =  arrayfun(@(z) onda_incidente(z), zz) ;
u_tot = u + u_inc;



figure()
subplot(3,3,1)
imagesc(x,x,real(u_inc))
title("R(u_{inc})")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(V, '-or')
% plot(Nodes, '+r', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')

subplot(3,3,2)
imagesc(x,x,imag(u_inc))
title("I(u_{inc})")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')

subplot(3,3,3)
imagesc(x,x,abs(u_inc))
title("|u_{inc}|")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')



subplot(3,3,4)
imagesc(x,x,real(u))
title("R(u_{scat})")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')
% 

subplot(3,3,5)
imagesc(x,x,imag(u))
title("I(u_{scat})")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')



subplot(3,3,6)
imagesc(x,x,abs(u))
title("|u_{scat}|")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')


subplot(3,3,7)
imagesc(x,x,real(u_tot))
title("R(u_{tot})")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')


subplot(3,3,8)
imagesc(x,x,imag(u_tot))
title("I(u_{tot})")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')



subplot(3,3,9)
imagesc(x,x,abs(u_tot))
title("|u_{tot}|")
set(gca, 'YDir', 'normal')   % mantiene l'orientamento matematico
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V_collection)
% plot(Nodes, 'or', MarkerFaceColor='r')
% plot(Midpoints, '+b', MarkerFaceColor='b')


MyGiffer(xx,yy,u,V_collection,"u_scatt_sol_fond.gif")
MyGiffer(xx,yy,u_inc,V_collection, "u_inc_sol_fond.gif")
MyGiffer(xx,yy,u_tot, V_collection,"u_tot_sol_fond.gif")

