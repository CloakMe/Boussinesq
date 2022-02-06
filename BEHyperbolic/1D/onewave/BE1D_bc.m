clear;clc;
% constants
% st ; end
%-20 ; 15
startx=-40; endx = 40;
tau = 0.025;  h = 0.025;  sgm=1/2;  x(:,1) = startx:h:endx;  t_end=3;
beta1=1.5;   beta2=0.5;  alpha=3;   c=2;
beta=beta1/beta2;
estep = 4;
sw=0; %switch on//off the grid shift (for 1 wave)

% necessary stuff
%gg=inline('bt*al*vv*vv + (bt -1)*vv','vv','al','bt');

sx = size(x,1);
dh=deltah(7);
dhs=dh*dh/h^2;
Idh=(h^2*eye(7)-dh);
    
    A = h^2*eye(7)/tau^2 - dh/tau^2 - sgm*dh + sgm*dhs;
    B = -2*A - dh + dhs;
    
    %B = -2*h^2*eye(sx)/tau^2 + (2/tau^2)*dh;
    a11=A(1,1);
    sA = A(3,1:5);
    sB = B(3,1:5);
    sdh = dh(2,1:3);
    sdhs = dhs(3,1:5);
    sIdh = Idh(2,1:3);
    clear A;clear B; clear dhs; clear dh; clear Idh;
    %1Bu CJLou ______________
    shift = 5;
    for l=1:sx
        v(l,1)= u_ex(x(l)+shift,0,c);% + u_ex(x(l)-shift,0,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % v(l,1)= u_ex(x(l)+shift,0,c);
    end
    %v(1,1)=0;v(end,1)=0;
    
    %2Pu CJLou ______________
    
    for l=1:sx  %g(u) u (du/dt)(x,0):
%        ggv0(l)=gg(v(l,1),alpha,beta);
        v1(l) =dudt_ex(x(l)+shift,0,c);% + dudt_ex(x(l)-shift,0,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % v1(l) =dudt_ex(x(l)+shift,0,c);
    end

    %b=(dh*v(:,1)) - (dhs*v(:,1)) + (dh*ggv0');
    %y=Idh\b; % /h^2 se sykrashtawa
    
    for i=1:sx % v Ha 2Puq CJLou
        %v(i,2) = v(i,1) + tau*v1(i) + tau^2*y(i)/2;
        v(i,2) = u_ex(x(i,1)+shift,tau,c);
    end
            
    figure(1)
    plot(x,v1,'g',x,v(:,1),'b',x,v(:,2),'r');%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         
    t(1)=0;t(2)=tau;
    k=2;
    %yoyo = tau
    s=1;e=1;
    E(1)=0;II(1)=0;

    vmo=v(2:end-1,1); vz = v(2:end-1,2);
    woc1 = 1;
    woc2 = 0;
    tic
    tau2=2*tau;
    while(t(k)<t_end)
        t(k+1) = t(k) + tau;
        bndry(1,1) = u_ex(shift+startx-h,t(k-1),c); bndry(1,2) = u_ex(shift+startx,t(k-1),c); bndry(1,3) = u_ex(shift+endx,t(k-1),c); bndry(1,4) = u_ex(shift+endx+h,t(k-1),c); 
        bndry(2,1) = u_ex(shift+startx-h,t(k),c);   bndry(2,2) = u_ex(shift+startx,t(k),c);   bndry(2,3) = u_ex(shift+endx,t(k),c);   bndry(2,4) = u_ex(shift+endx+h,t(k),c); 
        bndry(3,1) = u_ex(shift+startx-h,t(k+1),c); bndry(3,2) = u_ex(shift+startx,t(k+1),c); bndry(3,3) = u_ex(shift+endx,t(k+1),c); bndry(3,4) = u_ex(shift+endx+h,t(k+1),c); 
        
        ldtn(1) = sA(1)*bndry(3,1) + sA(2)*bndry(3,2);
        ldtn(2) =                    sA(1)*bndry(3,2);
        ldtn(3) =                                    sA(5)*bndry(3,3);
        ldtn(4) =                                    sA(4)*bndry(3,3) + sA(5)*bndry(3,4);
    
        Avmo(1) = sA(1)*bndry(1,1) + sA(2)*bndry(1,2);
        Avmo(2) =                    sA(1)*bndry(1,2);
        Avmo(3) =                                    sA(5)*bndry(1,3);
        Avmo(4) =                                    sA(4)*bndry(1,3) + sA(5)*bndry(1,4);
        
        Bvz(1) = sB(1)*bndry(2,1) + sB(2)*bndry(2,2);
        Bvz(2) =                    sB(1)*bndry(2,2);
        Bvz(3) =                                     sB(5)*bndry(2,3);
        Bvz(4) =                                     sB(4)*bndry(2,3) + sB(5)*bndry(2,4);
        %A*vpo == dh*g1 - A*vmo - B*vz
        for ii=1:4  
            pbad(ii) = Avmo(ii) + Bvz(ii);  
        end
        
        csol = SIT_bc2(sA,sB,sdh,bndry,pbad,ldtn,vz,vmo,sx,alpha,beta);
        yoyo=t(k)
        v(1,k+1) = bndry(3,2);
        v(2:end-1,k+1) = csol;
        v(end,k+1) = bndry(3,3);
        if(mod(k,estep)==0)
            II(e)=sum(v(:,k))*h;
            EL(e)=LE(v(:,k-1),v(:,k),v(:,k+1),sdh,sIdh,h,tau,sgm);
            E(e) = EL(e) + NLE(v(:,k),v(:,k+1),sx,alpha,beta);
            e=e+1;
        end
       
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
            
            k=k+1;
             vmo = v(2:end-1,k-1); 
             vz  = v(2:end-1,k); 
        else
           k=k+1;  
           
            vmo = vz;
            vz  = v(2:end-1,k);
        end
    end
    %E(1) = E(2);
    toc
    figure(2)
    for l=1:sx
         uex2(l) =  u_ex(x(l)+shift,t_end,c);% + u_ex(x(l)-shift,t_end,-1.5);%<^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        %uex2(l) =  u_ex(x(l)+shift,t_end,c);
    end
    plot(x(:,end),v(:,1),'c',x(:,end),v(:,2),'g',x(:,end),v(:,end),'r',x(:,end),uex2,'k');
    if(size(E,2)*size(E,1)>6)
        EE = E - floor(E(1));
        figure(3)
        %plot(t(2:estep:end-1),E,'g',t(2:estep:end-1),E2,'b');
        plot(t(2:estep:end-1),EE,'c');
        title('Energy');
        dE = max(E(1:end)) - min(E(1:end))
        E(end)
        figure(5)
        plot(t(2:estep:end-1),II,'c');
        title('Integral');
    end
       
    sol_size = size(v)
   
 
    maxvalv = 1.05*max(v(:,end));
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
        Tr = u_ex(x(:,ff)+shift,t(j),c);% + u_ex(x(:,ff)-shift,t(j),-1.5);
        plot(x(:,ff),v(:,j),x(:,ff),Tr,'k');%+h*(ff-1)
        axis([startx x(end,end) (-0.06) maxvalv]);%maxvalv]);
        F(:,j) = getframe;
        %Error = max(abs(v(:,j)-Tr'))
     end
     movie(F,1)
     
     
     %load v_sol;
     %dlmwrite('vbc_t10_h025.dat', v)
     
     %vbc_t10_h025e=v(:,end);
     %save vbc_t10_h025.mat vbc_t10_h025e t x;
     Error = max(abs(v(:,end)-Tr'))
    