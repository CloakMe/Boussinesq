function [niv] = SIT_v11(sA,sB,sdh,a11,b11,vz,vmo,sx,al,bt)
%===============================================
    %vzl = [vz(2:end); 0];
    %vzr = [0; vz(1:end-1)];
    %g1 = bt*al*(vzl.^2 + 10*vz.^2 + vzr.^2)/12;
    
    g1 = bt*al*vz.^2;

     % A = h^2*eye(sx)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    pb = BMM(sA,vmo',a11) + BMM(sB,vz',b11);
    b = BMM(sdh,g1') - pb;
    
   % y=L\b;
    niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa

    