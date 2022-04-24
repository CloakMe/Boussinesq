function [niv] = SIT_v2(Idh,Itaudh,dh,dhh,AuxM,vz,vmo,sx,alpha,beta)
    eps = 10^(-13);
    gg1=inline('bt*al*(vvu^2 + vvu*vvd +vvd^2)/3 + (bt-1)*(vvu+vvd)/2','vvd','vvu','al','bt');
    %spparms('bandden',0.5);
%===============================================
    for l=1:sx
        g1(l)=gg1(vmo(l),vz(l),alpha,beta);
    end
   % y=L\b;
    w = tridiagsolve(Idh,dh*g1'); % /h^2 se sykrashtawa
    
    pb = AuxM*vz + Itaudh*vmo;
    b = w - pb;
    niv0 = tridiagsolve(Itaudh,b);
    %w0(1:20)
    for l=1:sx
       g1(l)=gg1(vmo(l),niv0(l),alpha,beta);
    end
    %y=L\b;
    w = tridiagsolve(Idh,dh*g1'); % /h^2 se sykrashtawa
    
    pb = AuxM*vz + Itaudh*vmo;
    b = w - pb;
    niv = tridiagsolve(Itaudh,b);
%===============================================
    
   % max(abs(w0-w))
    mmax=max(abs(niv));
    cnt=0;
    while(max(abs(niv0 - niv))>eps*mmax)
        niv0=niv;
        mmax=max(abs(niv));
        if(cnt>15)
           fprintf(1,'More than 15 it...')
           pause
        end
        cnt=cnt+1; %7
%===============================================
        for l=1:sx
            g1(l)=gg1(vmo(l),niv(l),alpha,beta);
        end
        %y=L\b;
        w = tridiagsolve(Idh,dh*g1'); % /h^2 se sykrashtawa
    
        pb = AuxM*vz + Itaudh*vmo;
        b = w - pb;
        niv = tridiagsolve(Itaudh,b);
        %max(A*w + b - dh*g1')
%===============================================
    end
    if(cnt < 6 && 8 < cnt)
        cnt
    end


    
    
    