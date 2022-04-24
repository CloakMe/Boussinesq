function [niv] = SIT_vsgm(sA,sB,sdh,a11,b11,vz,vmo,h,al,bt,sgm)
    eps = 10^(-13);
%    gg1=inline('bt*al*(vvu^2 + vvu*vvd +vvd^2)/3 + (bt-1)*(vvu+vvd)/2','vvd','vvu','al','bt');
    %spparms('bandden',0.5);
    
    g1 = (bt*al)*((1/12) * vmo.^2 + (11/12)* vz.^2 ) + (bt-1)*( (1/12)*vmo + (11/12)*vz );
    g1 = vz - BMM(sdh,vz')'/h^2 + g1;
    
    b =  BMM(sdh,g1');
    
   % y=L\b;
    niv0 = 2*vz - vmo +  pentsolv(a11,sA,b); % /h^2 se sykrashtawa
    
    g1 = (bt*al)*((1/12) * vmo.^2 + (10/12)* vz.^2 + (1/12) * niv0.^2 ) + (bt-1)*( (1/12)*vmo + (10/12)*vz + (1/12) * niv0 );
    g1 = vz - BMM(sdh,vz')'/h^2 + g1;
    
    b =  BMM(sdh,g1');
    %y=L\b;
    niv = 2*vz - vmo +  pentsolv(a11,sA,b); % /h^2 se sykrashtawa
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
        g1 = (bt*al)*((1/12) * vmo.^2 + (10/12)* vz.^2 + (1/12) * niv0.^2 ) + (bt-1)*( (1/12)*vmo + (10/12)*vz + (1/12) * niv0 );
        g1 = vz - BMM(sdh,vz')'/h^2 + g1;
    
        b =  BMM(sdh,g1');
    
        %y=L\b;
        niv = 2*vz - vmo + pentsolv(a11,sA,b); % /h^2 se sykrashtawa
        %max(A*niv + b - dh*g1')
%===============================================
        %max(abs(niv0-niv));
    end



