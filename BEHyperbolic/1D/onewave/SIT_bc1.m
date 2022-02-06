function [niv] = SIT_bc1(sA,sB,sdh,bndry,pbad,ldtn,vz,vmo,sx,alpha,beta,tt)
    startx=-20; endx = 15; c=2;
    hh = 0.05;  sgm=1/2;  xx = startx-hh:hh:endx+hh;
    shift = 5;

    a11 = sA(3);
    eps = 10^(-13);
    gg1=inline('bt*al*(vvu^2 + vvu*vvd +vvd^2)/3 + (bt-1)*(vvu+vvd)/2','vvd','vvu','al','bt');
    for l=1:size(xx,2)
        %g1(l)=gg1(vmo(l),vz(l),alpha,beta);
        u_tochpo(l)=u_ex(shift+xx(l),tt+0.05,c);
        u_tochz(l)=u_ex(shift+xx(l),tt      ,c);
        u_tochmo(l)=u_ex(shift+xx(l),tt-0.05,c);        
    end
    %spparms('bandden',0.5);
    sx=sx-2;
    sx
    size(vmo)
    size(vz)
%===============================================
    g1ad(1) = gg1(bndry(1,2),bndry(3,2),alpha,beta);
    Tg1ad(1) = gg1(u_ex(shift+xx(2),tt-0.05,c),u_ex(shift+xx(2),tt+0.05,c),alpha,beta);
    for l=1:sx
        %g1(l)=gg1(vmo(l),vz(l),alpha,beta);
        g1(l)=gg1(vmo(l),u_ex(shift+xx(l+2),tt+0.05,c),alpha,beta); %u_ex(shift+xx(l+2),tt+0.5,c)
    end
    g1ad(2) = gg1(bndry(1,3),bndry(3,3),alpha,beta);
    Tg1ad(2) = gg1(u_ex(shift+xx(end-1),tt-0.05,c),u_ex(shift+xx(end-1),tt+0.05,c),alpha,beta);
    norm((g1ad - Tg1ad),2)
    for l=1:sx+2
        %g1(l)=gg1(vmo(l),vz(l),alpha,beta);
        Tg1(l)=gg1(u_tochmo(l+1),u_tochpo(l+1),alpha,beta);
    end
   
       %A*vpo == dh*g1 - A*vmo - B*vz
    pb(1) = sA(3:5)*vmo(1:3) + sB(3:5)*vz(1:3) + pbad(1); 
    b(1) = sdh(1)*g1ad(1) + sdh(2:3)*g1(1:2)'  - pb(1);
    pb(2) = sA(2:5)*vmo(1:4) + sB(2:5)*vz(1:4) + pbad(2);
    b(2) = sdh*g1(1:3)' - pb(2);
    for i=1:sx-4
        pb(i+2) = sA*vmo(i:i+4) + sB*vz(i:i+4);
        b(i+2) = sdh*g1((i+1):(i+3))' - pb(i+2);
    end
    pb(sx-1) = sA(1:4)*vmo(sx-3:sx) + sB(1:4)*vz((sx-3):sx) + pbad(3);
    b(sx-1) = sdh*g1((sx-2):sx)' - pb(sx-1);
    pb(sx)   = sA(1:3)*vmo((sx-2):sx) + sB(1:3)*vz((sx-2):sx) + pbad(4);
    b(sx)   = sdh(1:2)*g1((sx-1):sx)' + sdh(3)*g1ad(2) - pb(sx);
    sTG1 = size(Tg1)
    size(xx,2)
    size(u_tochmo)
    size(u_tochz)
    for i=1:(size(xx,2)-4)
        Tb(i) = sdh*Tg1((i):(i+2))' - sA*u_tochmo(i:i+4)' - sB*u_tochz(i:i+4)';
    end
    norm_b = norm(b-Tb,2)
    figure(6)
    plot(xx(3:end-2),b-Tb);
    figure(7)
    plot(xx(3:end-2),g1-Tg1(2:end-1));
    figure(8)
    plot(xx(3:end-2),vz'-u_tochz(3:end-2));
    
    b(1) = b(1) - ldtn(1);    b(2) = b(2) - ldtn(2);
    b(sx-1) = b(sx-1) - ldtn(3)- ldtn(3);    b(sx) = b(sx) - ldtn(4);
    
   % y=L\b;
    niv0 = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    figure(5)
    plot(xx(3:end-2),niv0,'g',xx(2:end-1),u_ex(5+xx(2:end-1),tt+0.05,c));
    
    %niv0(1:20)
    g1ad(1) = gg1(bndry(1,2),bndry(3,2),alpha,beta);
    for l=1:sx
       g1(l)=gg1(vmo(l),niv0(l),alpha,beta);
    end
    g1ad(2) = gg1(bndry(1,3),bndry(3,3),alpha,beta);
    
    b(1) = sdh(1)*g1ad(1) + sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
    for i=1:sx-2
        b(i+1) = sdh*g1(i:i+2)' - pb(i+1);
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
        %g1ad(1) = gg1(bndry(1,2),bndry(3,2),alpha,beta);
        for l=1:sx
            g1(l)=gg1(vmo(l),niv(l),alpha,beta);
        end
        %g1ad(2) = gg1(bndry(1,3),bndry(3,3),alpha,beta);
        
        b(1) = sdh(1)*g1ad(1) + sdh(2)*g1(1) + sdh(3)*g1(2) - pb(1);
        for i=1:sx-2
            b(i+1) = sdh*g1(i:i+2)' - pb(i+1);
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



