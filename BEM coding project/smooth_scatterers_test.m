close all
clear
%% scatterer liscio e confronto nel caso del cerchio

a = -1;
b = 1;
n_points = 150;


k = 30;
x = linspace(a,b,n_points);
[xx,yy] = meshgrid(x,x);                
zz = xx + 1i * yy; 
L = 50;  % serie in Z la approssimiamo come serie da -L a L
phi = pi/6;
rr = abs(zz);
theta = angle(zz); % in radianti

l = -L:L;
R = 0.25;
J = arrayfun(@(x) besselj(x,k*R), l, 'UniformOutput',false);

H =  arrayfun(@(x) besselh(x, k*rr)./besselh(x, k*R), l, 'UniformOutput',false);


u_scat_multipole = zeros(size(zz));
for j = 1:length(l)
    a_l = 1i.^l(j) * exp(-1i*l(j)*phi)*cell2mat(J(j));

    u_scat_multipole = u_scat_multipole - a_l.*cell2mat(H(j)) .*exp(1i*l(j)*theta);
end

% mettiamo maschera perchè ''dentro'' il disco la soluzione esplode
mask_outside = rr > R;
u_scat_multipole(~mask_outside) = 0;


% calcolo errore come differenza tra u_scat multipoli e u_scat collocazione
ins= 10:10:300;
errori = zeros(size(ins));
for mm = 1:length(ins)
    nr_punti = ins(mm);
    t = linspace(0, 2*pi, nr_punti);

    
    % circonferenza come in 4.21
    x = cos(t)*(0.25);
    y = sin(t)*(0.25);
    
    
    p = polyshape(x(1:end-1),y(1:end-1), 'SolidBoundaryOrientation','ccw');
    plot(p, "FaceColor",[1,1,1], "EdgeColor", [1,1,1], "FaceAlpha",1)
    hold on
    axis equal
    
    
    verts = p.Vertices;
    
    V = {verts(:,1) + 1i*verts(:,2)}; % e ora lo passiamo al nostro solutore solito
    
   
    h_ideale = pi/(5*k); %così evito di non beccare le oscillazioni che sono legate a k (remark 5.9) ;
    
    q_quadratura = ceil(3*k); % idem
    
    [Nodes, Midpoints, Hj, Dirs] = create_nodes(V, h_ideale);
    
    d = [sqrt(3)/2, 0.5];
    
    
    plane_wave = @(z) cos(k*(real(z)*d(1)+imag(z)*d(2))) + 1i*sin(k*(real(z)*d(1)+imag(z)*d(2))); % exp(1i*k*(real(z)*d(1) + imag(z)*d(2)));
    x0 = 1.5+1.5i;
    sol_fondamentale_x0 = @(z) 0.25i*besselh(0,k*abs(z-x0));
    
    
    onda_incidente = plane_wave;
    
    u =build_2D_solution(V, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    
    
    u_inc =  arrayfun(@(z) onda_incidente(z), zz) ;
    u_tot = u + u_inc;
    
    u(~mask_outside) = 0;
    dx = xx(1,2) - xx(1,1);
    dy = yy(2,1)-yy(1,1);
    num = sqrt(sum(sum(abs(u - u_scat_multipole).^2 ))* dx * dy);
    den = sqrt(sum(sum(abs(u_scat_multipole).^2)) * dx * dy);
    errori(mm) = num./den;
    disp("fatto con m" + string(mm))
end
   

% grafico per test multipoli
u_tot_multipole = u_inc + u_scat_multipole;
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
plot_polygon(V)

subplot(3,3,2)
imagesc(x,x,imag(u_inc))
title("I(u_{inc})")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)

subplot(3,3,3)
imagesc(x,x,abs(u_inc))
title("|u_{inc}|")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)


subplot(3,3,4)
imagesc(x,x,real(u_scat_multipole))
title("R(u_{scat})")
set(gca, 'YDir', 'normal')   
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)
subplot(3,3,5)
imagesc(x,x,imag(u_scat_multipole))
title("I(u_{scat})")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)



subplot(3,3,6)
imagesc(x,x,abs(u_scat_multipole))
title("|u_{scat}|")
set(gca, 'YDir', 'normal') 
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)


subplot(3,3,7)
imagesc(x,x,real(u_tot_multipole))
title("R(u_{tot})")
set(gca, 'YDir', 'normal')   
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)


subplot(3,3,8)
imagesc(x,x,imag(u_tot_multipole))
title("I(u_{tot})")
set(gca, 'YDir', 'normal')  
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)


subplot(3,3,9)
imagesc(x,x,abs(u_tot_multipole))
title("|u_{tot}|")
set(gca, 'YDir', 'normal') 
colormap("parula")
colorbar
axis equal
hold on
plot_polygon(V)


 
exportgraphics(fig, "smooth_scatterer_multipole.png", 'Resolution', 300);

% per le animazioni
% % % MyGiffer(xx,yy,u,V,"u_scatt_smooth.gif")
% % % MyGiffer(xx,yy,u_inc,V, "u_inc_smooth.gif")
% % % MyGiffer(xx,yy,u_tot, V,"u_tot_smooth.gif")
% % 





