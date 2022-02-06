clear;clc;
% constants
% st ; end
%-20 ; 15
startx=-20;
tau = 0.05;  h = 0.05;  sgm=1/2;  x(:,1) = startx:h:20;  t_end=1;
beta1=1.5;   beta2=0.5;  alpha=3;
beta=beta1/beta2;
sw=0; %switch on//off the grid shift (for 1 wave)

% necessary stuff
gg=inline('bt*al*vv*vv + (bt -1)*vv','vv','al','bt');

sx = size(x,1);
dh=deltah_v2(sx);
dhs=dh*dh/h^2;
Idh=(h^2*eye(sx)-dh);
    %1Bu CJLou ______________
    
    for l=1:sx
        v(l,1)= u_ex(x(l)+5,0,2);% + u_ex(x(l)-5,0,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % v(l,1)= u_ex(x(l)+5,0,2);
    end
    %v(1,1)=0;v(end,1)=0;
    
    %2Pu CJLou ______________
    
    for l=1:sx  %g(u) u (du/dt)(x,0):
        ggv0(l)=gg(v(l,1),alpha,beta);
        v1(l) =dudt_ex(x(l)+5,0,2);% + dudt_ex(x(l)-5,0,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % v1(l) =dudt_ex(x(l)+5,0,2);
    end

    b=(dh*v(:,1)) - (dhs*v(:,1))+ (dh*ggv0');
    y=Idh\b; % /h^2 se sykrashtawa
    
    for i=1:sx % v Ha 2Puq CJLou
        v(i,2) = v(i,1) + tau*v1(i) + tau^2*y(i)/2;
    end
    figure(1)
    plot(x,v1,'g',x,v(:,1),'b',x,v(:,2),'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    A = h^2*eye(sx)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    B = -2*A - dh + dhs;
    [L,U]=lu(A);
    t(1)=0;t(2)=tau;
    k=2;
    yoyo = tau
    s=1;
    vmo=v(:,1); vz = v(:,2);
    woc1 = 1;
    woc2 = 0;
    tic
    tau2=2*tau;
    while(t(k)<t_end)
        %v(:,k+1) = SIT_v4(sA,sB,sdh,a11,vz,vmo,sx,alpha,beta);
        v(:,k+1) = SIT_lu1(A,B,L,U,dh,vz,vmo,sx,alpha,beta);
        %vt = (v(:,k+1) - v(:,k))/tau;
        %vec1=(-dh/h^2)\vt;
        %EL(k)=(vec1')*vt + (vt')*vt + tau^2*(sgm-1/4)*(((Idh*vt)')*vt)/h^2 +...  %
        %    (( v(:,k) + v(:,k+1) + (-dh)*(v(:,k)+ v(:,k+1))/h^2)')*(v(:,k) + v(:,k+1))/4;
        %E(k) = EL(k) + NLE(v(:,k),v(:,k+1),sx,alpha,beta);
       
        if(mod(t(k),woc1)==woc2 && sw == 1)
            woc1 = 1;woc2 = 0;
            s=s+1;
            [a,b] = max(v(:,k+1));
            midd=floor(sx/2);
            shift=b-midd;
            %plot(x(:,1),v(:,k+1),'k',x(:,1),v(:,k),'g')
            temp = v(:,k+1);
            v(1:sx-2*shift,k+1)=temp(2*shift+1:sx);
            v(sx-2*shift+1:sx,k+1)=v(sx-2*shift,k+1)*v(sx-2*shift+1:sx,2)/v(sx-2*shift+1,2);
            
            temp = v(:,k);
            v(1:sx-2*shift,k)=temp(2*shift+1:sx);
            v(sx-2*shift:sx,k)=v(sx-2*shift-1,k)*v(sx-2*shift:sx,2)/v(sx-2*shift+1,2);
            
            x(:,s) = x(:,s-1) + h*(2*shift);
            xshift(s) = h*(2*shift);
            
            %figure(2);
            %plot(x(:,1),v(:,k+1),'k',x(:,1),v(:,k),'g')
            k=k+1;
             t(k)=(k-1)*tau;
             yoyo=t(k)
             vmo = v(:,k-1); 
             vz  = v(:,k); 
            % max(abs(vmo - vz))
        else
           k=k+1;  
           
            t(k)=(k-1)*tau;
            yoyo=t(k)
            %vmt = vmo;
            vmo = vz;
            vz  = v(:,k);
           % max(abs(vmo - vz))
        end
    end
    %E(1) = E(2);
    toc
    figure(2)
    for l=1:sx
         uex2(l) =  u_ex(x(l)+5,t_end,2);% + u_ex(x(l)-5,t_end,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        %uex2(l) =  u_ex(x(l)+5,t_end,2);
    end
    plot(x(:,end),v(:,1),'c',x(:,end),v(:,2),'g',x(:,end),v(:,end),'r',x(:,end),uex2,'k');
   % figure(3)
   % plot(t(1:end-1),E,'g');
   % title('Energy');
       
    xiks_size = size(x)
    time_size = size(t)
    sol_size = size(v)
   
 
maxvalv = 1.5*max(v(:,end));
if(isnan(maxvalv))
    maxvalv = 1.5*max(v(:,1));
    if(isnan(maxvalv))
        maxvalv = 5;
    end
end

ff=1;
woc1 = 1; woc2 = 0;
figure(4)
     for j = 1:size(t,2) 
        if(mod(t(j),woc1) == woc2 && t(j)~=0 && t(j)~=t_end && sw==1)
            woc1 = 1;woc2 = 0;
            ff=ff+1;
        end
        Tr = u_ex(x(:,ff)+5,t(j),2);% + u_ex(x(:,ff)-5,t(j),-1.5);
        plot(x(:,ff),v(:,j),x(:,ff),Tr,'k');%+h*(ff-1)
        axis([startx x(end,end) (-0.3) maxvalv]);
        F(:,j) = getframe;
        %Error = max(abs(v(:,j)-Tr'))
     end
     movie(F,1)
     
     
     %load v_sol;
     dlmwrite('v_t10_h0125', v)
     
     v_t10_h0125_end=v(:,end);
     save v_t10_h0125.mat t x v_t10_h0125_end;
     Error = max(abs(v(:,end)-Tr'))
    