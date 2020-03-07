function [psi,kapa]= error_rate(u1,u2,u3)
% we find the error psi and the order of convergence kapa using Runge's method
% here u1, u2 and u3 are the solution over embeded grids

Y2=Comp(u2);
psi1=max(abs(u1-Y2));

Y3=Comp(u3);
psi2=max(abs(u2-Y3));

kapa=log2(psi1/psi2);
psi=psi1*psi1/(psi1-psi2);

end