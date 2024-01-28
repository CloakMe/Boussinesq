function [vu,dvu,v,tt,II] = BE1D_tv3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0,hOrder,tauOrder)
%A solution using Taylor Series Approach
%The function returns the solution "v", the time derivative "dtv" and the times "tt" for which the
%solution "v" has been saved
% dE = max(E) - min(E);
% dII = max(II) - min(II);
% Here E stands for the energy values and II stands for the integral values along
% different time layers "tt".
if(nargin == 12)
    hOrder = 2;
end
if(nargin == 13)
    tauOrder = 4;
end
beta = beta1/beta2;
x = start_x:h:end_x;
sx = size(x,2);
    if(hOrder == 2)
        dh = deltah(7);
        sdh = dh(2,1:3);
        sdh11 = dh(2,2);
    else
        dh = deltaOh4(7);
        sdh = dh(3,1:5);
        sdh11 = dh(3,3);
    end
    Idh = h^2*eye(7) - deltaOh4(7);
    sIdh = Idh(3,1:5);
    sIdh11 = Idh(1,1);
   
    v1 = u_t0';
    dv1  = dudt_t0';
    [d2v1, d3v1, d4v1, d5v1] = calc_der_sub(v1,dv1,sdh,sIdh,sdh11,sIdh11,h,alpha,beta);
        
    v2  = v1  + tau*dv1  + tau^2*d2v1/2 + (tauOrder == 4) * tau^3*d3v1/6 + (tauOrder == 4) * tau^4*d4v1/24; 
    dv2 = dv1 + tau*d2v1 + tau^2*d3v1/2 + (tauOrder == 4) * tau^3*d4v1/6 + (tauOrder == 4) * tau^4*d5v1/24;
    
    %figure(1)
    %plot(x,dv1,'g',x,dv2,'y',x,v1,'b',x,v2,'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    %figure(1)
    %plot(x,v1,'b',x,dv1,'g',x,d2v1,'y',x,d3v1,'r',x,d4v1,'k');
    %xlabel('x');    ylabel('u');
    t(1)=0;t(2)=tau;
    E(1)=0;dE=0;dII = 0;
    k=2;
    s=1; e=1;
    vmo=v1; vz = v2;
    dvmo=dv1; dvz = dv2;
    tt = zeros(1,t_end*10);
    II = zeros(1,t_end*10);
    v = zeros(sx,10*t_end-1);
    dtv = zeros(sx,10*t_end-1);

    while(t(k)<t_end)
        
        [d2vz, d3vz, d4vz, d5vz] = calc_der_sub(vz,dvz,sdh,sIdh,sdh11,sIdh11,h,alpha,beta);
    
        vu  =  vz +  tau*dvz + tau^2*d2vz/2 + (tauOrder == 4) * tau^3*d3vz/6 + (tauOrder == 4) * tau^4*d4vz/24; 
        dvu = dvz + tau*d2vz + tau^2*d3vz/2 + (tauOrder == 4) * tau^3*d4vz/6 + (tauOrder == 4) * tau^4*d5vz/24;
        
        %vu = SIT_m1(invdh,dh,tau,vz,vmo,alpha,beta);

        if(mod(k,estep)==0)
            tt(e)=k*tau;
            v(:,e) = vz;
            dtv(:,e) = dvz;
            II(e)=sum(vz)*h;
            %EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            %E(e) = EL(e) + NLE(vz,vu,sx,h,alpha,beta);
            if(abs(tt(e)-t_end) <= 2*tau)
                TEND = tt(e)
            end
            e=e+1;
        end
        if( abs( ( t_end-tau + 1.00e-010) - t(k))  < tau )
             if(abs(tt(end)-t_end) > 1.00e-010)
                warning('t_end missed!'); 
                tt(end)
                t(end)      
             else
                 SOL_SAVED = 1
             end
        end
          k=k+1;  
        if(max(abs(vu))>50)
            warning('solution diverges!'); 
           return; 
        end
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k))) < 0.9*tau)
                yoyo=t(k)
            end
            vmo = vz;     dvmo = dvz;
            vz  = vu;     dvz  = dvu;
    end

    dII = norm(II - II(end),2);
    %dE = norm(E - E(end),2);
        %[d2vz, d3vz, d4vz, d5vz] = calc_der(vz,dvz,sdh,sIdh,h,sx,alpha,beta);
        %time= - 0.013189885977779;
        %vta_t10_h025_t025e  =  vz +  time*dvz + time^2*d2vz/2 + time^3*d3vz/6 + time^4*d4vz/24;

       % save vta_t10_h025_t025.mat vta_t10_h025_t025e tt x;
end
    
 

    