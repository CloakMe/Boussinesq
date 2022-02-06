function [niv] = SIT_bc2(sA,sB,sdh,bndry,pbad,ldtn,vz,vmo,sx,al,bt)
    a11 = sA(3);
    eps = 10^(-13);
    gg1=inline('bt*al*(vvu^2 + vvu*vvd +vvd^2)/3 + (bt-1)*(vvu+vvd)/2','vvd','vvu','al','bt');
    %startx=-20; endx = 10; c=2;
    %hh = 0.1;  sgm=1/2;  xx = startx+hh:hh:endx-hh;

    %spparms('bandden',0.5);
    sx=sx-2;
 %===============================================
    %g1ad(1) = gg1(bndry(1,2),bndry(2,2),al,bt);
    g1ad(1) = bt*al*(bndry(2,2)^2 + bndry(2,2)*bndry(1,2) +bndry(1,2)^2)/3 + (bt-1)*(bndry(2,2)+bndry(1,2))/2;
    g1 = bt*al*(vz.^2 + vz'*vmo +vmo.^2)/3 + (bt-1)*(vz+vmo)/2;
    %g1ad(2) = gg1(bndry(1,3),bndry(2,3),al,bt);
    g1ad(2) = bt*al*(bndry(2,3)^2 + bndry(2,3)*bndry(1,3) +bndry(1,3)^2)/3 + (bt-1)*(bndry(2,3)+bndry(1,3))/2;
    
       %A*vpo == dh*g1 - A*vmo - B*vz
       %pb = A*vmo + B*vz;
       %b = dh*g1' - pb;
    pb(1) = sA(3:5)*vmo(1:3) + sB(3:5)*vz(1:3) + pbad(1); 
    b(1) = sdh(1)*g1ad(1) + sdh(2:3)*g1(1:2)  - pb(1);
    pb(2) = sA(2:5)*vmo(1:4) + sB(2:5)*vz(1:4) + pbad(2);
    b(2) = sdh*g1(1:3) - pb(2);
    for i=1:sx-4
        pb(i+2) = sA*vmo(i:i+4) + sB*vz(i:i+4);
        b(i+2) = sdh*g1((i+1):(i+3)) - pb(i+2);
    end
    pb(sx-1) = sA(1:4)*vmo(sx-3:sx) + sB(1:4)*vz((sx-3):sx) + pbad(3);
    b(sx-1) = sdh*g1((sx-2):sx) - pb(sx-1);
    pb(sx)   = sA(1:3)*vmo((sx-2):sx) + sB(1:3)*vz((sx-2):sx) + pbad(4);
    b(sx)   = sdh(1:2)*g1((sx-1):sx) + sdh(3)*g1ad(2) - pb(sx);
    
    b(1) = b(1) - ldtn(1);    b(2) = b(2) - ldtn(2);
    b(sx-1) = b(sx-1) - ldtn(3);    b(sx) = b(sx) - ldtn(4);
    
   % y=L\b;
    niv0 = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    %niv0(1:20)
    g1ad(1) = bt*al*(bndry(1,2)^2 + bndry(1,2)*bndry(1,2) +bndry(1,2)^2)/3 + (bt-1)*(bndry(1,2)+bndry(1,2))/2;
    g1 = bt*al*(niv0.^2 + niv0.*vmo +vmo.^2)/3 + (bt-1)*(niv0+vmo)/2;
    g1ad(2) = bt*al*(bndry(3,3)^2 + bndry(3,3)*bndry(1,3) + bndry(1,3)^2)/3 + (bt-1)*(bndry(3,3)+bndry(1,3))/2;
    
    b(1) = sdh(1)*g1ad(1) + sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
    for i=1:sx-2
        b(i+1) = sdh*g1(i:i+2) - pb(i+1);
    end
    b(sx) = sdh(1)*g1(sx-1) + sdh(2)*g1(sx) + sdh(3)*g1ad(2) - pb(sx);
    
    b(1) = b(1) - ldtn(1);    b(2) = b(2) - ldtn(2);
    b(sx-1) = b(sx-1) - ldtn(3);    b(sx) = b(sx) - ldtn(4);

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
        %g1ad(1) = bt*al*(bndry(1,2)^2 + bndry(1,2)*bndry(1,2) +bndry(1,2)^2)/3 + (bt-1)*(bndry(1,2)+bndry(1,2))/2;
        g1 = bt*al*(niv0.^2 + niv0.*vmo +vmo.^2)/3 + (bt-1)*(niv0+vmo)/2;
        %g1ad(2) = bt*al*(bndry(3,3)^2 + bndry(3,3)*bndry(1,3) + bndry(1,3)^2)/3 + (bt-1)*(bndry(3,3)+bndry(1,3))/2;
        
        b(1) = sdh(1)*g1ad(1) + sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
        for i=1:sx-2
            b(i+1) = sdh*g1(i:i+2) - pb(i+1);
        end
        b(sx) = sdh(1)*g1(sx-1) + sdh(2)*g1(sx) + sdh(3)*g1ad(2) - pb(sx);
        
        b(1) = b(1) - ldtn(1);    b(2) = b(2) - ldtn(2);
        b(sx-1) = b(sx-1) - ldtn(3);    b(sx) = b(sx) - ldtn(4);
    
        %y=L\b;
        niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
     %   niv(1:20)
    end



