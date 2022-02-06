function [v,tt,dE,dII] = BE1D_vn1(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0)
%A solution using Numerov type of numerical scheme!
%The function returns the solution "v" and the times "tt" for which the
%solution "v" has been saved
% dE = max(E) - min(E);
% dII = max(II) - min(II);
% Here E stands for the energy values and II stands for the integral values along
% different time layers "tt".

x = start_x:h:end_x; 
beta = beta1/beta2;
sgm  = 0;
sx = size(x,2);
dh=deltah(7);
dhs=dh*dh/h^2;
Idh=(h^2*eye(7)-dh);
    %tldaIdh = h^2*eye(sx)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    tldaIdh = h^2*eye(7)/tau^2 - dh/tau^2;
    dhIdh = dh*(h^2*eye(7)-dh)/h^2;
    
    A = tldaIdh - dhIdh/12;
    B = 2*tldaIdh + (5/6)*dhIdh;

    a11=A(1,1);
    sA = A(3,1:5);
    b11 = B(1,1);
    sB = B(3,1:5);
    sdh = dh(2,1:3);
    sdhs = dhs(3,1:5);
    dhs11= dhs(1,1);
    sIdh = Idh(2,1:3);
    sIdh11 = Idh(1,1);
    sdh11 = dh(1,1);
    %clear A;clear B; clear dhs; clear dh; clear Idh;
    

    v1 = u_t0';
    dv1 = dudt_t0';
    ggv0 = beta*alpha*v1.^2 + (beta-1)*v1;
    
    b(1) = (sdh(2:3)*v1(1:2)) - dhs11*v1(1) - sdhs(4:5)*v1(2:3)+ (sdh(2:3)*ggv0(1:2));
    b(2) = (sdh*v1(1:3)) - (sdhs(2:5)*v1(1:4))+ (sdh*ggv0(1:3));
    for l=3:sx-2
        b(l)=(sdh*v1(l-1:l+1)) - (sdhs*v1(l-2:l+2))+ (sdh*ggv0(l-1:l+1));
    end
    b(sx-1) = (sdh*v1(sx-2:sx)) - (sdhs(1:4)*v1(sx-3:sx))+ (sdh*ggv0(sx-2:sx));
    b(sx) = sdh(1:2)*v1(sx-1:sx) - sdhs(1:2)*v1(sx-2:sx-1) - dhs11*v1(sx)+ (sdh(1:2)*ggv0(sx-1:sx));
    y=diag3solv(sIdh,b); % /h^2 se sykrashtawa
    
    
    v2 = v1 + tau*dv1 + tau^2*y/2;
    
    figure(1)
    plot(x,dv1,'g',x,v1,'b',x,v2,'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
     
    t(1)=0;t(2)=tau;
    E(1)=0;II(1) = 0;
    k=2;
    s=1;e=1;
    vmo=v1; vz = v2;

    tic

    while(t(k)<t_end)
        %vu = SIT_vn4(sA,sB,sdh,a11,b11,vz,vmo,sx,alpha,beta);
        vu = SIT_vn1(sA,sB,sdh,a11,b11,vz,vmo,sx,alpha,beta);
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            v(:,e) = vu;
            II(e)=sum(vz)*h;
            EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            E(e) = EL(e) + NLE(vz,vu,sx,alpha,beta);
            e=e+1;
        end
 
           k=k+1;  
           
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k)))<=5*tau/6)
                yoyo=t(k)
            end
            vmo = vz;
            vz  = vu;
    end
    
    toc

    if(size(E,2)*size(E,1)>6)
        EE = E - floor(E(1));
        figure(3)
        %plot(t(2:estep:end-1),E,'g',t(2:estep:end-1),E2,'b');
        plot(t(2:estep:end-1),E,'c');
        title('Energy');
        dE = max(E) - min(E);
        figure(5)
        plot(t(2:estep:end-1),II,'c');
        title('Integral');
        dII = max(II) - min(II);
    end

     
    