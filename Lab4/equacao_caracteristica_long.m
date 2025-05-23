clear all 
syms Xu Xw Xq w0 g theta_0 Zu Zw Zq u0 phi_0 Mu Mw Mq

A_long_syms = [Xu, Xw, Xq-w0, -g*cos(theta_0);
          Zu, Zw, Zq+u0, -g*cos(phi_0)*sin(theta_0);
          Mu, Mw, Mq, 0;
          0, 0, cos(phi_0), 0];
      
coeficients = charpoly(A_long_syms);

A = coeficients(1)
B = coeficients(2)
C = coeficients(3)
D = coeficients(4)
E = coeficients(5)