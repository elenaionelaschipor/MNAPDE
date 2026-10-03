function f = sol_fondamentale(x,y,k)
f = 0.25i * besselh(0, 1, k*abs(x-y));
end