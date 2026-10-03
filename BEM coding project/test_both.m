clear
close all


a = -1;
b = 2;
n_points = 200;
V_collection = {[0, 1, 1i]};

k = 20;
h_ideale = pi/(10*k); %così evito di non beccare le oscillazioni che sono legate a k (remark 5.9) ;
q_quadratura = ceil(3*k); % idem


% Ns = find_ideal_mesh(V_collection, h_ideale);
[Nodes, Midpoints, Hj, Dirs] = create_nodes(V_collection, h_ideale);

d = [0.5,sqrt(3)/2];


plane_wave = @(z) cos(k*(real(z)*d(1)+imag(z)*d(2))) + 1i*sin(k*(real(z)*d(1)+imag(z)*d(2))); % exp(1i*k*(real(z)*d(1) + imag(z)*d(2)));
x0 = 1.5+1.5i;
sol_fondamentale_x0 = @(z) 0.25i*besselh(0,k*abs(z-x0));


onda_incidente = plane_wave;


tempi_galerkin_triangolo = zeros(10,1);
tempi_collocazione_triangolo = tempi_galerkin_triangolo;


disp("le triangle")
for i = 1:length(tempi_galerkin_triangolo)
    tic;
    u_gal = build_2D_solution_galerkin(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    tempi_galerkin_triangolo(i) = toc;
    
    tic;
    u_col = build_2D_solution(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    tempi_collocazione_triangolo(i)=toc;
    
    disp("fatto i = " + string(i))

end


tempi_galerkin_cerchio = zeros(10,1);
tempi_collocazione_cerchio= tempi_galerkin_cerchio;

disp("cerchioooooooo")

t = linspace(0, 2*pi, 200); % Finer sampling gives a smoother edge
% circonferenza come in 4.21
x = cos(t)*(0.25);
y = sin(t)*(0.25);


p = polyshape(x(1:end-1),y(1:end-1), 'SolidBoundaryOrientation','ccw');
plot(p, "FaceColor",[1,1,1], "EdgeColor", [1,1,1], "FaceAlpha",1)
hold on
axis equal

verts = p.Vertices;

V_collection = {verts(:,1) + 1i*verts(:,2)}; % e ora lo passiamo al nostro solutore solito!!!


h_ideale = pi/(5*k); %così evito di non beccare le oscillazioni che sono legate a k (remark 5.9) ;

q_quadratura = ceil(3*k); % idem

for i = 1:length(tempi_galerkin_cerchio)
    tic;
    u_gal = build_2D_solution_galerkin(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    tempi_galerkin_cerchio(i) = toc;
    tic;
    u_col = build_2D_solution(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    tempi_collocazione_cerchio(i)=toc;
    disp("fatto i " + string(i))
end



disp(" 2 triangoli, cambio onda incidente ")
V_collection = {[0,1,1i],[2,3,3+1i]};
onda_incidente = sol_fondamentale_x0;



tempi_galerkin_2tri = zeros(10,1);
tempi_collocazione_2tri= tempi_galerkin_2tri;

for i = 1:length(tempi_galerkin_2tri)
    tic;
    u_gal = build_2D_solution_galerkin(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    tempi_galerkin_2tri(i) = toc;
    tic;
    u_col = build_2D_solution(V_collection, h_ideale,k, onda_incidente,a,b,n_points, q_quadratura);
    tempi_collocazione_2tri(i)=toc;

    disp("Fatto i = " + string(i))
end