function [niv] = SIT_vn12(sA,sB,sdh,a11,b11,vz,vmo,al,bt)
    eps = 10^(-13);
    %gg=inline('bt*al*vv^2 + (bt-1)*vv','vv','al','bt');
    %ggn = inline('vu/12 + vn*5/6 + vd/12','vu','vn','vd');
    %spparms('bandden',0.5);
%===============================================
    
    ggz = bt*al*vz.^2 + (bt-1)*vz;
    ggmo = bt*al*vmo.^2 + (bt-1)*vmo;
    g1 = ggz/12 + (5/6)*ggz + ggmo/12;
   
    pb = BMM(sA,vmo',a11) + BMM(sB,vz',b11);
    b = BMM(sdh,g1') - pb;
    
    %y=L\b;
    niv0 = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    %niv0(1:20)
    ggz = bt*al*vz.^2 + (bt-1)*vz;
    ggmo = bt*al*vmo.^2 + (bt-1)*vmo;
    ggpo = bt*al*niv0.^2 + (bt-1)*niv0;
    g1 = ggpo/12 + ggz*5/6 + ggmo/12;
    
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
        ggz = bt*al*vz.^2 + (bt-1)*vz;
        ggmo = bt*al*vmo.^2 + (bt-1)*vmo;
        ggpo = bt*al*niv0.^2 + (bt-1)*niv0;
        g1 = ggpo/12 + ggz*5/6 + ggmo/12;
        
        b = BMM(sdh,g1') - pb;
        
        %y=L\b;
        niv = pentsolv(a11,sA,b); % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
     %   niv(1:20)
    end

