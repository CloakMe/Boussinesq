function [niv] = SIT_v4(sA,sB,sdh,a11,b11,vz,vmo,sx,al,bt)
    eps = 10^(-13);
%    gg1=inline('bt*al*(vvu^2 + vvu*vvd +vvd^2)/3 + (bt-1)*(vvu+vvd)/2','vvd','vvu','al','bt');
    %spparms('bandden',0.5);
%===============================================
    g1 = bt*al*(vmo.^2 + vmo.*vz +vz.^2)/3 + (bt-1)*(vmo+vz)/2;

     % A = h^2*eye(sx)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    pb = BMM(sA,vmo',a11) + BMM(sB,vz',b11);
    b = BMM(sdh,g1') - pb;
    
   % y=L\b;
    niv0 = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    
    g1 = bt*al*(niv0.^2 + niv0.*vmo + vmo.^2)/3 + (bt-1)*(niv0+vmo)/2;
    
    b = BMM(sdh,g1') - pb;
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
        g1 = bt*al*(niv.^2 + niv.*vmo + vmo.^2)/3 + (bt-1)*(niv+vmo)/2;
        b = BMM(sdh,g1') - pb;
    
        %y=L\b;
        niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
    end



