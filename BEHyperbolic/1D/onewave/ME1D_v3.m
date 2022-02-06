function [v,tt,dE,dII] = ME1D_v3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0)
%A solution using five diagonal matrices
% The approximation of the nonlinear term has the following presentation (IM):
% g_1=inline('(bt*al/3)*(v_up^2 + v_up*v_down +v_down^2) + ((bt-1)/2)*(v_up+v_down)','v_up','v_down','al','bt');
% The approximation inlcudes the unknown layer v_up which is found
% through some iterative procedure.
%The function returns the solution "v" and the times "tt" for which the
%solution "v" has been saved
% dE = max(E) - min(E);
% dII = max(II) - min(II);
% Here E stands for the energy values and II stands for the integral values along
% different time layers "tt".
beta = beta1/beta2;
x = start_x:h:end_x;
sx = size(x,2);
dh = deltah(sx)/h^2;

    new_sdh = h^2*dh(2,1:3);
    new_sdh11 = h^2*dh(1,1);
    Idh = h^2*eye(7) - deltah(7);
    sIdh = Idh(2,1:3);
    sIdh11 = Idh(1,1);
    invdh = inv(inv(dh) - eye(sx));
    
    v1 = u_t0';
    
    dv1  = dudt_t0';
    d2v1 = (invdh*(alpha*beta*v1.*v1 + (beta-1)*v1) + dh*v1);
    d3v1 = (invdh*(2*alpha*beta*dv1.*v1 + (beta-1)*dv1) + dh*dv1);
    d4v1 = (invdh*(2*alpha*beta*(dv1.*dv1 + v1.*d2v1) + (beta-1)*d2v1) + dh*d2v1);
    d5v1 = (invdh*(2*alpha*beta*(3*dv1.*d2v1 + v1.*d3v1) + (beta-1)*d3v1) + dh*d3v1);
    
    v2 = v1 + tau*dv1 + tau^2*d2v1/2 + tau^3*d3v1/6 + tau^4*d4v1/24; 
    dv2 = dv1 + tau*d2v1 + tau^2*d3v1/2 + tau^3*d4v1/6 + tau^4*d5v1/24;
    
    figure(1)
    plot(x,dv1,'g',x,dv2,'y',x,v1,'b',x,v2,'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    
    t(1)=0;t(2)=tau;
    E(1)=0;II(1) = 0;
    k=2;
    s=1; e=1;
    vmo=v1; vz = v2;
    dvz = dv2; %dvmo=dv1;
    v = zeros(sx,20*t_end-1);
    %dtv = zeros(sx,20*t_end-1);
    
    tic
    while(t(k)<t_end)
        
        d2vz = (invdh*(alpha*beta*vz.*vz + (beta-1)*vz) + dh*vz);
        d3vz = (invdh*(2*alpha*beta*dvz.*vz + (beta-1)*dvz) + dh*dvz);
        d4vz = (invdh*(2*alpha*beta*(dvz.*dvz + vz.*d2vz) + (beta-1)*d2vz) + dh*d2vz);
        d5vz = (invdh*(2*alpha*beta*(3*dvz.*d2vz + v1.*d3vz) + (beta-1)*d3vz) + dh*d3vz);
    
        vu = vz + tau*dvz + tau^2*d2vz/2 + tau^3*d3vz/6 + tau^4*d4vz/24; 
        dvu = dvz + tau*d2vz + tau^2*d3vz/2 + tau^3*d4vz/6 + tau^4*d5vz/24;
        
        %vu = SIT_m1(invdh,dh,tau,vz,vmo,alpha,beta);

        if(mod(k,estep)==0)
            tt(e)=k*tau;
            v(:,e) = vu;
            %dtv(:,e) = dvu;
            II(e)=sum(vz)*h;
            EL(e)=LE(vmo,vz,vu,new_sdh,sIdh,new_sdh11,sIdh11,h,tau,sgm);
            E(e) = EL(e) + NLE(vz,vu,sx,alpha,beta);
            e=e+1;
        end
          k=k+1;  
           
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k)))==0)
                yoyo=t(k)
            end
            vmo = vz;     %dvmo = dvz;
            vz  = vu;     dvz  = dvu;
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
    dE=0;dII=0;
   
    