clear all 
syms Yv Yp w0 g phi_0 Yr u0 Lv_prime Lp_prime Lr_prime theta_0 Nv_prime Np_prime Nr_prime

A_lat_syms = [Yv, Yp + w0, g*cos(phi_0)*cos(phi_0), Yr - u0;
        Lv_prime, Lp_prime, 0, Lr_prime;
        0, 1, 0, cos(phi_0)*tan(theta_0);
        Nv_prime, Np_prime, 0, Nr_prime];
      
coeficients = charpoly(A_lat_syms);

A = coeficients(1)
B = coeficients(2)
C = coeficients(3)
D = coeficients(4)
E = coeficients(5)