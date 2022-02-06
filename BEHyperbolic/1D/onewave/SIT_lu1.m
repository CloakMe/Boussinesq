function [niv] = SIT_lu1(A,B,L,U,dh,vz,vmo,sx,alpha,beta)
    eps = 10^(-13);
    gg1=inline('bt*al*(vvu^2 + vvu*vvd +vvd^2)/3 + (bt-1)*(vvu+vvd)/2','vvd','vvu','al','bt');
    %spparms('bandden',0.5);
%===============================================
    for l=1:sx
        g1(l)=gg1(vmo(l),vz(l),alpha,beta);
    end
        
    pb = A*vmo + B*vz;
    b = dh*g1' - pb;
    y=L\b;
    niv0 = U\b; % /h^2 se sykrashtawa
    %niv0(1:20)
    for l=1:sx
       g1(l)=gg1(vmo(l),niv0(l),alpha,beta);
    end
    b = dh*g1' - pb;
    y=L\b;
    niv = U\b; % /h^2 se sykrashtawa
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
        for l=1:sx
            g1(l)=gg1(vmo(l),niv(l),alpha,beta);
        end
        b = dh*g1' - pb;
        y=L\b;
        niv = U\b; % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
     %   niv(1:20)
    end

