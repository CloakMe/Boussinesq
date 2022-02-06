function [niv] = SIT_v6(sA,sB,sdh,a11,b11,vz,vmo,sx,al,bt)
    eps = 10^(-13);
    %G = inline('(al*bt*vv^3)/3 + ((bt-1)*vv^2)/2','vv','al','bt');
    %gg1=inline('2*(GF - GS)/(vvu-vvd)','GF','GS','vvu','vvd','al','bt');

%===============================================
     vapz = vz;
     GF = (al*bt*vapz.^3)/3 + ((bt-1)*vapz.^2)/2;
     vapz = (vmo+vz)/2;
     GS = (al*bt*vapz.^3)/3 + ((bt-1)*vapz.^2)/2;
     GFS = GF - GS;
        g1 = 2*GFS./(vz-vmo);

     % A = h^2*eye(sx)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    pb = BMM(sA,vmo',a11) + BMM(sB,vz',b11);
    b = BMM(sdh,g1') - pb;
        
   % y=L\b;
    niv0 = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    %niv0(1:20)
          
    vapz = (niv0+vz)/2;
    GF = (al*bt*vapz.^3)/3 + ((bt-1)*vapz.^2)/2;
    vapz = (vmo+vz)/2;
    GS = (al*bt*vapz.^3)/3 + ((bt-1)*vapz.^2)/2;
    GFS = GF - GS;
       g1 = 2*GFS./(niv0-vmo);
    
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
        if(cnt>25)
            error('ERROR; too much iterations in SIT! ');
        end
        cnt=cnt+1; %5
%===============================================
        vapz = (niv0+vz)/2;
        GF = (al*bt*vapz.^3)/3 + ((bt-1)*vapz.^2)/2;
        vapz = (vmo+vz)/2;
        GS = (al*bt*vapz.^3)/3 + ((bt-1)*vapz.^2)/2;
        GFS = GF - GS;
            g1 = 2*GFS./(niv0-vmo);

        b = BMM(sdh,g1') - pb;
    
        %y=L\b;
        niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
     %   niv(1:20)
    end



