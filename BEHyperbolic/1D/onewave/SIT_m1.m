function [niv] = SIT_m1(iA,dh,tau,vz,vmo,alpha,beta)
    eps = 10^(-13);
    %gg=inline('bt*al*vv^2 + (bt-1)*vv','vv','al','bt');
    %ggn = inline('vu/12 + vn*5/6 + vd/12','vu','vn','vd');
    %spparms('bandden',0.5);
%===============================================
    niv0 = vz;
    niv0 = 2*vz - vmo +  (tau^2/12)*( iA*(  alpha*beta*niv0.*niv0 + (beta-1)*niv0 +...
        10*(beta-1)*vz + alpha*beta*vz.*vz + 10*(beta-1)*vz +...
        (beta-1)*vmo + alpha*beta*vmo.*vmo + (beta-1)*vmo  ) + dh*(niv0 + 10*vz + vmo)  );

    %y=L\b;
    niv = 2*vz - vmo +  (tau^2/12)*( iA*(  alpha*beta*niv0.*niv0 + (beta-1)*niv0 +...
        10*(beta-1)*vz + alpha*beta*vz.*vz + 10*(beta-1)*vz +...
        (beta-1)*vmo + alpha*beta*vmo.*vmo + (beta-1)*vmo  ) + dh*(niv0 + 10*vz + vmo)  );
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
        niv = 2*vz - vmo +  (tau^2/12)*( iA*(  alpha*beta*niv0.*niv0 + (beta-1)*niv0 +...
                                    10*alpha*beta*vz.*vz + 10*(beta-1)*vz +...
                                    alpha*beta*vmo.*vmo + (beta-1)*vmo ) + dh*(niv0 + 10*vz + vmo)  );
        
%===============================================
        %max(abs(niv0-niv));
     %   niv(1:20)
    end

