function [dv2,dv3,dv4,dv5,dv6,dv7,dv8]=calc_der_rs(r,v,dv1,gamma1,gamma2,order)

if(order<4 && order>8)
    error('order value is wrong!');
end

dv2 = - v / r + gamma1 * v + gamma2 * v^2;

dv3 = - dv1 / r + v / r ^ 2 + gamma1 * dv1 + 2 * gamma2 * v * dv1;

dv4 = - dv2 / r + 2 * dv1 / r ^ 2 - 2 * v / r ^ 3 + gamma1 * dv2 + 2 * gamma2 * dv1 ^ 2 + 2 * gamma2 * v * dv2;
if(order == 4) dv5=0; dv6=0;  dv7=0; dv8=0; return; end
    
dv5 = - dv3 / r + 3 * dv2 / r ^ 2 - 6 * dv1 / r ^ 3 + 6 * v / r ^ 4 + gamma1 * dv3 + 6 * gamma2 * dv1 * dv2 + 2 * gamma2 * v * dv3;
if(order == 5) dv6=0;  dv7=0; dv8=0; return; end

dv6 = - dv4 / r + 4 * dv3 / r ^ 2 - 12 * dv2 / r ^ 3 + 24 * dv1 / r ^ 4 - 24 * v / r ^ 5 + gamma1 * dv4 + 6 * gamma2 * dv2 ^ 2 + 8 * gamma2 * dv1 * dv3 + 2 * gamma2 * v * dv4;
if(order == 6)  dv7=0; dv8=0; return; end

dv7 = - dv5 / r + 5 * dv4 / r ^ 2 - 20 * dv3 / r ^ 3 + 60 * dv2 / r ^ 4 - 120 * dv1 / r ^ 5 + 120 * v / r ^ 6 + gamma1 * dv5 + 20 * gamma2 * dv2 * dv3 + 10 * gamma2 * dv1 * dv4 + 2 * gamma2 * v * dv5;
if(order == 7)  dv8=0; return; end

dv8 = - dv6 / r + 6 * dv5 / r ^ 2 - 30 * dv4 / r ^ 3 + 120 * dv3 / r ^ 4 - 360 * dv2 / r ^ 5 + 720 * dv1 / r ^ 6 - 720 * v / r ^ 7 + gamma1 * dv6 + 20 * gamma2 * dv3 ^ 2 + 30 * gamma2 * dv2 * dv4 + 12 * gamma2 * dv1 * dv5 + 2 * gamma2 * v * dv6;
