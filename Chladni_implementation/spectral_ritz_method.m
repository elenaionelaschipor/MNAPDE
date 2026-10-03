% definiamo le u_m
clear
close all

s = 6;
N = s+1;
syms x y
k_vals = zeros(1, N+1);
for m = 0:N
    if mod(m,2)==0
        k_vals(m+1) = vpasolve(tan(x)+tanh(x)==0, x, m*pi/2 - pi/4);
    else
        k_vals(m+1) = vpasolve(tan(x)-tanh(x)==0, x, (m-1/2)*pi/2);
    end
end

function u_m = calcola_u(m,x, k_vals)

if mod(m,2) == 0
    % m pari       
    % soluzione (numerica) di f1 = 0      
    
    if m == 0
        u_m = 1/sqrt(2);
    else
        k_m = k_vals(m+1); 
        u_m = (cosh(k_m)*cos(k_m*x) + cos(k_m)*cosh(k_m*x))/(sqrt(cosh(k_m)^2 + cos(k_m)^2));
    end

elseif mod(m,2) == 1 
    if m ==1
        u_m = sqrt(3/2)*x;
    else
        k_m = k_vals(m+1); 
        u_m = (sinh(k_m)*sin(k_m*x) + sin(k_m)*sinh(k_m*x))/(sqrt(sinh(k_m)^2 - sin(k_m)^2));
    end
end
end


% costuiamo la matrice K 


K_tilde = zeros((N+1)^2,(N+1)^2);

i = 0;
mu = 0.225;
ux = sym(zeros(1,N+1)) ; %zeros(1,s+2);
dux = sym(zeros(1,N+1));
d2ux = sym(zeros(1,N+1));

uy = sym(zeros(1,N+1)) ; %zeros(1,s+2);
duy = sym(zeros(1,N+1));
d2uy = sym(zeros(1,N+1));
   
for m=0:N
    ux(m+1) = calcola_u(m,x,k_vals);
    dux(m+1) = diff(ux(m+1), x);
    d2ux(m+1) = diff(dux(m+1), x);
    uy(m+1) = subs(ux(m+1),x,y);
    duy(m+1) = subs(dux(m+1), x,y);
    d2uy(m+1) =subs(d2ux(m+1),x,y);
end

for m=0:N
    for n=0:N
        i = i+1;
        j = 0;
        for p = 0:N
            for q = 0:N
                j = j+1;
                if mod(m+p, 2) == 0 && mod(n+q,2) == 0
                    u_p_x = ux(p+1);
                    d_u_p_x = dux(p+1);
                    d2_u_p_x = d2ux(p+1);

                    u_m_x = ux(m+1);
                    d_u_m_x = dux(m+1);
                    d2_u_m_x= d2ux(m+1);

                    u_n_y = uy(n+1);
                    d_u_n_y = duy(n+1);
                    d2_u_n_y = d2uy(n+1);

                    u_q_y = uy(q+1);
                    d_u_q_y = duy(q+1);
                    d2_u_q_y = d2uy(q+1);
                    

                    disp("sto lavorando con m,n,p,q" + string(m) + string(n) + string(p) + string(q))
                   
                    integranda = d2_u_m_x*u_n_y*d2_u_p_x*u_q_y ...
                                +u_m_x*d2_u_n_y*u_p_x*d2_u_q_y ...
                                + 2*mu*d2_u_m_x*u_n_y*u_p_x*d2_u_q_y ...
                                + 2*(1-mu)*d_u_m_x*d_u_n_y*d_u_p_x*d_u_q_y;
                    

                    integranda_num = matlabFunction(integranda, "Vars",{x,y});
                    
                    K_tilde(i,j) = integral2(integranda_num, -1,1,-1,1, "Vectorized",false);
                    %K_tilde(i,j) = int(int( integranda, x,-1,1), y, -1,1);

                end
            end
        end
    end
end

K = (K_tilde+K_tilde')/2;


[eigenvectors,D] = eig(K);
eigenvalues = diag(D);

eigenfunctions = sym(zeros(size(eigenvalues)));
% 
% for l = 1:N      % per ogni autovettore
%     w_l = sym(0);
%     for m = 0:s
%         for n = 0:s
%             k = m*(s+1) + n + 1;   % indice lineare per (m,n)
%             w_l = w_l + eigenvectors(k,l)*ux(m+1)*uy(n+1);
%         end
%     end
%     eigenfunctions(l) = w_l;
% end

for i = 1:(N+1)^2
    j = 0;
    eigenfunc = sym(0);
    for m =0:N
        for n=0:N
            j=j+1;
            eigenfunc = eigenfunc + eigenvectors(j,i)*ux(m+1)*uy(n+1);
        end
    end
    eigenfunctions(i) = eigenfunc;
end


% 
% close all
% for l = 1:N
%     figure()
%     fcontour(eigenfunctions(l), [-1 1 -1 1])
%     axis square 
% end