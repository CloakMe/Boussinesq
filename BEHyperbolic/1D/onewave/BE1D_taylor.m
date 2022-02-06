function [v,dtv,tt,E,II] = BE1D_taylor(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0)
%A solution using Taylor Series Approach
%The function returns the solution "v", the time derivative "dtv" and the times "tt" for which the
%solution "v" has been saved
% dE = max(E) - min(E);
% dII = max(II) - min(II);
% Here E stands for the energy values and II stands for the integral values along
% different time layers "tt".
beta = beta1/beta2;
x = start_x:h:end_x;
sx = size(x,2);

    dh = deltah(7);
    sdh = dh(2,1:3);
    Idh = h^2*eye(7) - deltah(7);
    sIdh = Idh(2,1:3);
    sIdh11 = Idh(1,1);sdh11 = dh(1,1);
    v1 = u_t0';
    dv1  = dudt_t0';
    [d2v1, d3v1, d4v1, d5v1] = calc_der(v1,dv1,sdh,sIdh,sdh11,sIdh11,h,alpha,beta);
        
    v2  = v1  + tau*dv1  + tau^2*d2v1/2 + tau^3*d3v1/6 + tau^4*d4v1/24; 
    dv2 = dv1 + tau*d2v1 + tau^2*d3v1/2 + tau^3*d4v1/6 + tau^4*d5v1/24;
    
    figure(1)
    plot(x,dv1,'g',x,dv2,'y',x,v1,'b',x,v2,'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    
    t(1)=0;t(2)=tau;
    E(1)=0;II(1) = 0;
    k=2;
    s=1; e=1;
    vmo=v1; vz = v2;
    dvmo=dv1; dvz = dv2;
    v = zeros(sx,20*t_end-1);
    dtv = zeros(sx,20*t_end-1);
    
    tic
    while(t(k)<t_end)
        
        [d2vz, d3vz, d4vz, d5vz] = calc_der(vz,dvz,sdh,sIdh,sdh11,sIdh11,h,alpha,beta);
    
        vu  =  vz +  tau*dvz + tau^2*d2vz/2 + tau^3*d3vz/6 + tau^4*d4vz/24; 
        dvu = dvz + tau*d2vz + tau^2*d3vz/2 + tau^3*d4vz/6 + tau^4*d5vz/24;
        
        %vu = SIT_m1(invdh,dh,tau,vz,vmo,alpha,beta);

        if(mod(k,estep)==0)
            tt(e)=k*tau;
            v(:,e) = vu;
            dtv(:,e) = dvu;
            II(e)=sum(vz)*h;
            EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            E(e) = EL(e) + NLE(vz,vu,sx,h,alpha,beta);
            if(abs(tt(e)-t_end) <= 2*tau)
                TEND = tt(e)
            end
            e=e+1;
        end
          k=k+1;  
           
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k))) < 0.9*tau)
                yoyo=t(k)
            end
            vmo = vz;     dvmo = dvz;
            vz  = vu;     dvz  = dvu;
    end
    toc
        %[d2vz, d3vz, d4vz, d5vz] = calc_der(vz,dvz,sdh,sIdh,h,sx,alpha,beta);
        %time= - 0.013189885977779;
        %vta_t10_h025_t025e  =  vz +  time*dvz + time^2*d2vz/2 + time^3*d3vz/6 + time^4*d4vz/24;

       % save vta_t10_h025_t025.mat vta_t10_h025_t025e tt x;
    %{
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
 %}
   
    