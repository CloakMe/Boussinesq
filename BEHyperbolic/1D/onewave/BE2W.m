clear;clc;
% constants
startx=-20;
tau = 0.1;  h = 0.1;  sgm=1/2;  x= startx:h:15;  t_end=2.5;
beta1=1.5;   beta2=0.5;  alpha=3;
beta=beta1/beta2;

% necessary stuff
gg=inline('bt*al*vv*vv + (bt -1)*vv','vv','al','bt');

sx = size(x,2);
dh=deltah(sx);
dhs=dh*dh/h^2;

    %1Bu CJLou ______________
    
    for l=1:sx
        v(l,1)= u_ex(x(l)+5,0,2,alpha,beta1,beta2);% + u_ex(x(l)-5,0,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % v(l,1)= u_ex(x(l)+5,0,2);
    end
    %v(1,1)=0;v(end,1)=0;
    
    %2Pu CJLou ______________
    
    for l=1:sx  %g(u) u (du/dt)(x,0):
        ggv0(l)=gg(v(l,1),alpha,beta);
        v1(l) =dudt_ex(x(l)+5,0,2,alpha,beta1,beta2);% + dudt_ex(x(l)-5,0,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % v1(l) =dudt_ex(x(l)+5,0,2);
    end

    b=(dh*v(:,1)) - (dhs*v(:,1)) + (dh*ggv0');
    AA=(h^2*eye(sx)-dh);
    y=AA\b; % /h^2 se sykrashtawa
    
    for i=1:sx % v Ha 2Puq CJLou
        v(i,2) = v(i,1) + tau*v1(i) + tau^2*y(i)/2;
    end

    for i=1:(sx)
       % vmt(i) = v(i,1) + (-tau)*v1(i) + tau^2*y(i)/2;
    end

    plot(x,v1,'g',x,v(:,1),'b');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    A = h^2*eye(sx)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    B = -2*A - dh + dhs;
    %[L,U]=lu(A);
    %Q=chol(-dh/h^2);
    %vec1=(Q*v(:,1));
    %EL(1)=(vec1')*vec1 + (v(:,1)')*v(:,1) + tau^2*(sgm-1/4)*(((AA*v(:,1))')*v(:,1))/h^2 +...
    %    (( v(:,1) + v(:,2) + (-dh)*(v(:,1)+ v(:,2))/h^2)')*(v(:,1) + v(:,2))/4;
    %E(1) = EL(1) + NLE(v(:,1),v(:,2),sx,alpha,beta);
    t(1)=0;t(2)=tau;
    k=2;
    vmo=v(:,1); vz = v(:,2);
    tic
    while(t(k)<t_end)
        v(:,k+1) = SIT(A,B,dh,vz,vmo,sx,alpha,beta);
        
        %vec1=(Q*v(:,k));
        %EL(k)=(vec1')*vec1 + (v(:,k)')*v(:,k) + tau^2*(sgm-1/4)*(((AA*v(:,k))')*v(:,k))/h^2 +...
        %    (( v(:,k) + v(:,k+1) + (-dh)*(v(:,k)+ v(:,k+1))/h^2)')*(v(:,k) + v(:,k+1))/4;
        %E(k) = EL(k) + NLE(v(:,k),v(:,k+1),sx,alpha,beta);
        
           k=k+1;  
           
        t(k)=(k-1)*tau;
        yoyo=t(k)
        %vmt = vmo;
        vmo = vz;
        vz  = v(:,k);
        max(abs(vmo - vz))
    end
    figure(1)
    for l=1:sx
         uex2(l) =  u_ex(x(l)+5,t_end,2);% + u_ex(x(l)-5,t_end,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        %uex2(l) =  u_ex(x(l)+5,t_end,2);
    end
    plot(x,v(:,1),'c',x,v(:,2),'g',x,v(:,end-2),'r',x,uex2,'k');
    %max(v(:,end))
    
   % figure(2)
   % plot(x,v(:,1),'c',x,v(:,end),'r',x,uex2,'k');
        
    xiks_size = size(x)
    time_size = size(t)
    sol_size = size(v)
   
    %dlmwrite('v_sol1', v)
    %load v_sol;
    %v=v_sol;
    figure(3)
    mesh(x',t',v')
    xlabel('x');    ylabel('t');    zlabel('u');
    %figure(7)
    %plot(t(1:end-1),E,'g');
    %title('Energy');
    %figure(5)
    %for l=1:sx
    %    uex2(l) =u_ex(x(l)+5,3.5,2) + u_ex(x(l)-5,3.5,-1.5);
    %end
   % plot(x,uex2,'k');
toc
figure(4)
     for j = 1:size(t,2) 
        Tr = u_ex(x+5,t(j),2);% + u_ex(x-5,t(j),-1.5);
        plot(x,v(:,j),x,Tr,'k');
        axis([-15 80 -0.5 2.6]);
        F(:,j) = getframe;
        j=j+1;
     end
     movie(F,1)
     
     %vv01=v(:,end);
     %save vv01.mat vv01 x t;
     Error = max(abs(v(:,end)-Tr'))
    