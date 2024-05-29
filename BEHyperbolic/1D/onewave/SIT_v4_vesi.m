function [niv, cnt] = SIT_v4_vesi(A,sigma,u1,u0,h,tau,alpha)
    eps = 10^(-13);
    %gg=inline('bt*al*vv^2 + (bt-1)*vv','vv','al','bt');
    %ggn = inline('vu/12 + vn*5/6 + vd/12','vu','vn','vd');
    %spparms('bandden',0.5);
%===============================================
   rightside=( 2*u1 - u0 - 2*sigma*tau*tau*vesi_deltah(u1,h) + ...
       sigma*tau*tau*vesi_deltah(u0,h) + ...
       2*sigma*tau*tau*vesi_delta2h(u1,h) - ...
       sigma*tau*tau*vesi_delta2h(u0,h) + ...
       tau*tau*(vesi_deltah(u1,h) - ...
       vesi_delta2h(u1,h) + ...
       alpha*vesi_deltah( (u0.^2 + u0.*u1 + u1.^2)/3, h)) )';
   
    niv0 = (A\rightside)';
    niv = niv0+eps*1000;

%     g1 = bt*al*(vz.^2 + vz.*vmo + vmo.^2)/3 + (bt-1)*(vz + vmo)/2;
%     g1 = vz - BMM(sdh,vz')'/h^2 + g1;
%     
%     b = BMM(sdh,g1');
%    
%     %y=L\b;
%     niv0 = 2*vz - vmo + pentsolv(a11,sA,b); % /h^2 se sykrashtawa
%     %niv0(1:20)
%     g1 = bt*al*(niv0.^2 + niv0.*vmo + vmo.^2)/3 + (bt-1)*(niv0 + vmo)/2;
%     g1 = vz - BMM(sdh,vz')'/h^2 + g1;
%     
%     b =  BMM(sdh,g1');
%     
%     %y=L\b;
%     niv = 2*vz - vmo + pentsolv(a11,sA,b); % /h^2 se sykrashtawa
%===============================================
    
   % max(abs(niv0-niv))
    mmax=max(abs(niv0));
    cnt=1;
    while(max(abs(niv0 - niv))>eps*mmax)
        niv0=niv;
        mmax=max(abs(niv));
        if(cnt>25)
            error('ERROR; too much iterations in SIT! ');
        end
        cnt=cnt+1; %5
        
        rightside=( 2*u1 - u0 - 2*sigma*tau*tau*vesi_deltah(u1,h) + ...
            sigma*tau*tau*vesi_deltah(u0,h) + ...
            2*sigma*tau*tau*vesi_delta2h(u1,h) - ...
            sigma*tau*tau*vesi_delta2h(u0,h) + ...
            tau*tau*(vesi_deltah(u1,h) - ...
            vesi_delta2h(u1,h) + ...
            alpha*vesi_deltah( (u0.^2 + u0.*niv0 + niv0.^2)/3, h)) )';

        niv = (A\rightside)';
%===============================================
%         g1 = bt*al*(niv0.^2 + niv0.*vmo + vmo.^2)/3 + (bt-1)*(niv0 + vmo)/2;
%         g1 = vz - BMM(sdh,vz')'/h^2 + g1;
%     
%         b =  BMM(sdh,g1');
%         
%         niv = 2*vz - vmo + pentsolv(a11,sA,b); % /h^2 se sykrashtawa
%===============================================
    end

    
    
    