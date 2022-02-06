function [vu,tt,v,E,II] = BE1D_v3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0)
%A solution using five diagonal matrices
% The approximation of the nonlinear term has the following presentation (IM):
% g_1=inline('(bt*al/3)*(v_up^2 + v_up*v_down +v_down^2) +
% ((bt-1)/2)*(v_up+v_down)','v_up','v_down','al','bt');
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
dh=deltah(7);
dhs=dh*dh/h^2;
Idh=(h^2*eye(7)-dh);
          
    v1 = u_t0';
    dv1 = dudt_t0';
    %------------------------------------------
    O4dh = deltaOh4(7);
    O4sdh = O4dh(3,1:5);
    O4sdh11 = O4dh(1,1);
    O4Idh = h^2*eye(7) - deltaOh4(7);
    O4sIdh = O4Idh(3,1:5);
    O4sIdh11 = O4Idh(1,1);

    %A = h^2*eye(7)/tau^2 - dh/tau^2;% - sgm*dh + sgm*dhs;
    A = h^2*eye(7)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    B = - dh + dhs;
    
    %B = -2*h^2*eye(sx)/tau^2 + (2/tau^2)*dh;
    a11=A(1,1);
    sA = A(3,1:5);
    b11 = B(1,1);
    sB = B(3,1:5);
    sdh = dh(2,1:3);
    sdh11 = dh(1,1);
    sdhs = dhs(3,1:5);
    dhs11= dhs(1,1);
    sIdh = Idh(2,1:3);
    sIdh11 = Idh(1,1);
    clear A;clear B; clear dhs; clear dh; clear Idh;
    
    [d2v1, d3v1, d4v1, d5v1] = calc_der(v1,dv1,O4sdh,O4sIdh,O4sdh11,O4sIdh11,h,alpha,beta);
        
    v2  = v1  + tau*dv1  + tau^2*d2v1/2 + tau^3*d3v1/6 + tau^4*d4v1/24; 
    
    %------------------------------------------------------
    
    %ggv0 = beta*alpha*v1.^2 + (beta-1)*v1;
    
    %b = BMM(sdh,v1') - BMM(sdhs,v1',dhs11) + BMM(sdh,ggv0');
    %y=diag3solv(sIdh,b); % /h^2 se sykrashtawa
        
    %v2 = v1 + tau*dv1 + tau^2*y/2;
    
    %figure(1)
    %plot(x,dv1,'g',x,v1,'b',x,v2,'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    t(1)=0;t(2)=tau;
    E(1)=0;II(1) = 0;dE=0;dII = 0;
    k=2;
    s=1; e=1;
    vmo=v1; vz = v2;
    v = zeros(sx,20*t_end-1);
    
    tic
    while(t(k)<t_end)
        %picardi g1:
        vu = SIT_v3(sA,sB,sdh,a11,b11,vz,vmo,h,alpha,beta);
        
        %Numerov approximation of g1:
        %vu = SIT_vn2(sA,sB,sdh,a11,b11,vz,vmo,sx,alpha,beta); 
        %No approximation of g1; no iterations, only using "v zero" inside g1:
        %vu = SIT_v10(sA,sB,sdh,a11,b11,vz,vmo,sx,alpha,beta);

        if(mod(k,estep)==0)
            tt(e)=k*tau;
            v(:,e) = vu;
            II(e)=sum(vz)*h;
            EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            E(e) = EL(e) + NLE(vz,vu,sx,h,alpha,beta);
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
           
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k)))==0)
                yoyo=t(k)
            end
            vmo = vz;
            vz  = vu;
    end
    toc

    