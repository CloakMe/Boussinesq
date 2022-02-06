function [niv] = SIT_vn4(sA,sB,sdh,a11,b11,vz,vmo,sx,al,bt)
    eps = 10^(-13);
    %gg=inline('bt*al*vv^2 + (bt-1)*vv','vv','al','bt');
    %ggn = inline('vu/12 + vn*5/6 + vd/12','vu','vn','vd');
    %spparms('bandden',0.5);
%===============================================
   
    ggz = bt*al*vz.^2 + (bt-1)*vz;
    ggmo = bt*al*vmo.^2 + (bt-1)*vmo;
    niv0 = vz;
    g1 = al*bt*(vz.^2 + (niv0.^2 - 2*niv0.*vmo + vmo.^2)/24 + vz.*(niv0 - 2*vz + vmo)/6) + (bt-1)*(niv0 + 10*vz + vmo)/12;
   
    pb(1) = (a11*vmo(1) + sA(4)*vmo(2) + sA(5)*vmo(3)) - (b11*vz(1) + sB(4)*vz(2) + sB(5)*vz(3));
    b(1) = sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
    pb(2) = sA(2:5)*vmo(1:4) - sB(2:5)*vz(1:4);
    b(2) = sdh*g1(1:3) - pb(2);
    for i=1:sx-4
        pb(i+2) = +sA*vmo(i:i+4) - sB*vz(i:i+4);
        b(i+2) = sdh*g1(i+1:i+3) - pb(i+2);
    end
    pb(sx-1) = sA(1:4)*vmo(end-3:end) - sB(1:4)*vz(end-3:end);
    b(sx-1) = sdh*g1(end-2:end) - pb(sx-1);
    pb(sx) = (sA(1)*vmo(end-2) + sA(2)*vmo(end-1) + a11*vmo(end)) - (sB(1)*vz(end-2) + sB(2)*vz(end-1) + b11*vz(end));
    b(sx) = sdh(1)*g1(end-1) + sdh(2)*g1(end) - pb(sx);
    
    %y=L\b;
    niv0 = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    %niv0(1:20)
    g1 = al*bt*(vz.^2 + (niv0.^2 - 2*niv0.*vmo + vmo.^2)/24 + vz.*(niv0 - 2*vz + vmo)/6) + (bt-1)*(niv0 + 10*vz + vmo)/12;

    
    b(1) = sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
    for i=1:sx-2
        b(i+1) = sdh*g1(i:i+2) - pb(i+1);
    end
    b(sx) = sdh(1)*g1(end-1) + sdh(2)*g1(end) - pb(sx);
    
    %y=L\b;
    niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
%===============================================
    
   % max(abs(niv0-niv))
    mmax=max(abs(niv0));
    cnt=0;
    while(max(abs(niv0 - niv))>eps*mmax)
        niv0=niv;
        mmax=max(abs(niv));
        if(cnt>15)
            error('ERROR; too much iterations in SIT! ');
        end
        cnt=cnt+1; %5
%===============================================
        g1 = al*bt*(vz.^2 + (niv0.^2 - 2*niv0.*vmo + vmo.^2)/24 + vz.*(niv0 - 2*vz + vmo)/6) + (bt-1)*(niv0 + 10*vz + vmo)/12;
 
        b(1) = sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
        for i=1:sx-2
            b(i+1) = sdh*g1(i:i+2) - pb(i+1);
        end
        b(sx) = sdh(1)*g1(end-1) + sdh(2)*g1(end) - pb(sx);
        
        %y=L\b;
        niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
     %   niv(1:20)
    end

