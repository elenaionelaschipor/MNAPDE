function P = SLP(z, Nodes, Dirs, Hj, s,w, k)


P = Hj.*sum(w.*besselh(0,1,k*abs(z - Nodes-Dirs.*Hj.*s)));



end