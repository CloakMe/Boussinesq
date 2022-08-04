function [tt, max_v, t, vpt]=BE2D_v5(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0)

delete F:\MATLAB6p5\work\MT\2D\SOL\*


beta = beta1/beta2;
sy = length(y);
sx = length(x);

max_v(1) = max(max(u_t0));
  
    dhb=deltah2d_Oh4(sx);
    [W,w] = eig(-dhb);
    
    DD = zeros(1,sx);
    IDH = DD;
    %I = eye(sx);

    for i = 1:sx
        %DD(i) = W(:,i)'*BMM([1 -4 1],W(:,i)')';%dhb*W(:,i);
        IDH(i) = W(:,i)'*BMM([(1/12) (-16/12) (h^2+5)  (-16/12) (1/12)],W(:,i)',(h^2+59/12))';%dhb*W(:,i);
    end
         
    vdah = zeros(sx,sy);
    b = vdah;
    VV = vdah;
    vz = vdah; vmo = vdah; vu = vdah;
    
    n_idh = 5;
    
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(u_t0,dudt_t0,W,IDH,n_idh,vdah,h,sx,al,beta,6);
    tau = tau/2;
    vmo = u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    dvmo = dudt_t0 + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vmo,dvmo,W,IDH,n_idh,vdah,h,sx,al,beta,6);
 
    vmo = vmo + tau*dvmo + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    dvmo = dvmo + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    
    
    
    
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vmo,dvmo,W,IDH,n_idh,vdah,h,sx,al,beta,6);
        
    vz = vmo + tau*dvmo + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    dvz = dvmo + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vz,dvz,W,IDH,n_idh,vdah,h,sx,al,beta,6);
        
    vz = vz + tau*dvz + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    dvz = dvz + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    
    
    
    
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vz,dvz,W,IDH,n_idh,vdah,h,sx,al,beta,6);
    
    vpo = vz + tau*dvz + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    dvpo = dvz + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vpo,dvpo,W,IDH,n_idh,vdah,h,sx,al,beta,6);
    
    vpo = vpo + tau*dvpo + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    %dvpo = dvpo + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    %====================================================
    tau = 2*tau;
    
    %[d2vz, d3vz, d4vz] = calc_der2d(vpo,dvpo,W,IDH,n_idh,vdah,h,sx,al,beta,4);
    
   % vu = vpo + tau*dvpo + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz ;
    %dvpo = dvpo + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    %figure(2)
    %mesh(x,y,(tau^2/2)*d2vz')
    %title('(tau^2/2)*d^2v/dt^2 , t=tau')
    %xlabel('x');            ylabel('y');
    
    %figure(2)
    %mesh(x,y,(O4dh_gmI(u_t0,0,vdah,vdah))/h^2');
    %title('DH * v')
    %xlabel('x');            ylabel('y');
    
    e=1;
    t(1)=0;t(2)=tau;
    k=2;
    vmt = u_t0;
    clear('u_t0');clear('dudt_t0');
    vdah = zeros(sx,sy);
    while(t(k)<t_end)

        %vzl = [vz(2:end); 0];
        %vzr = [0; vz(1:end-1)];
        %g1 = (bt*al*(vzl.^2 + 10*vz.^2 + vzr.^2) + (bt-1)*(vzl + 10*vz +  vzr))/12;
        vpt=beta*al*vz.*vz + (beta-1)*vz;
        b = W'*(O4dh_gmI(vpt,0,vdah)); %num_ap(vu,vdah)
    
        for j=1:sx
           VV(j,:) = pentsolv(IDH(j) - 1/12, [(1/12) (-16/12) IDH(j) (-16/12) 1/12],b(j,:));
        end
         
        vpt = -vmt + 16*vmo - 30*vz + 16*vpo - (12*tau^2/h^2)*(O4dh_gmI(vz,0,vdah)) - 12*tau^2*W*VV;
        %{
        figure(12);
        mesh(x,y,vpt');
        
        xlabel('x');            ylabel('y');
        %}
        IS_NaN = max(max(isnan(vpt)));
        max_v(k) = max(max(vpt));
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            %figure(2)
            %mesh(x,y,(O4dh_gmI(vz,0,vdah))');
            %title('DH * v')
            %xlabel('x');            ylabel('y');
            %figure(4)
            %mesh(x,y,((vpt - vmo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            
            str = ['SOL\MM_' num2str(tt(e)) '.mat'];
            save (str,  'vpt');%(1:xstep:end,1:xstep:end);
            %II(e)=sum(vz)*h;
            %EL(e)=LE(vmo,vz,vpt,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            %E(e) = EL(e) + NLE(vz,vpt,sx,alpha,beta);
            e=e+1;
        end
        if(IS_NaN > 0)
           warning('ERROR; BE2D NaN values; Stopping! ');
           return;
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
            vmt = vmo;
            vmo = vz;
            vz  = vpo;
            vpo = vpt;
    end
  
    
   clear('W');    clear('vmo'); 
    sol_size = size(vpt)

  
function vdah=dh_gmI(M,gm,vdah)   %dh operator if gm = 0

      j=1;
      vdah(j,:) = - (2+gm)*M(j,:) +  M(j+1,:);
      for j=2:size(M,1)-1
          vdah(j,:) =  M(j-1,:) - (2+gm)*M(j,:) +  M(j+1,:);
      end
      j=size(M,1);
      vdah(j,:) =  M(j-1,:) - (2+gm)*M(j,:) ;
      
      
      j=1;
      vdah(:,j) = vdah(:,j) - (2)*M(:,j) +  M(:,j+1);
      for j=2:size(M,1)-1
          vdah(:,j) = vdah(:,j)  +  M(:,j-1) - (2)*M(:,j) +M(:,j+1);
      end
      j=size(M,1);
      vdah(:,j) =  vdah(:,j) +M(:,j-1) - (2)*M(:,j) ;
      
function resmm=num_ap(M,vdah)
      vdah(1:end-1,:) = M(2:end,:); 
      vdah(2:end,:) = vdah(2:end,:) + M(1:end-1,:);
      resmm = (vdah + 10*M)/12;
      
function vdah=O4dh_gmI(M,gm,vdah)   %dh operator if gm = 0
    
      j=1;
      vdah(j,:) = - (2.5 -1/12 + gm)*M(j,:) + (16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      j=2;
      vdah(j,:) =    + (16/12)*M(j-1,:) - (2.5 + gm)*M(j,:) +  ( 16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      for j=3:size(M,1)-2
          vdah(j,:) = (-1/12)*M(j-2,:) + (16/12)*M(j-1,:) - (2.5 + gm)*M(j,:) +  (16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      end
      vdah(end-1,:) = (-1/12)*M(end-3,:)+ (16/12)*M(end-2,:) -  (2.5 + gm)*M(end-1,:) + (16/12)*M(end,:);
      vdah(end,:) =                        (-1/12)*M(end-2,:) + (16/12)*M(end-1,:) - (2.5 - 1/12 + gm)*M(end,:);
  
      
      j=1;
      vdah(:,j) = vdah(:,j) - (2.5 -1/12 )*M(:,j) + (16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      j=2;
      vdah(:,j) = vdah(:,j) +   (16/12)*M(:,j-1) - (2.5 )*M(:,j)  +  ( 16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      for j=3:size(M,2)-2
        vdah(:,j) = vdah(:,j) + (-1/12)*M(:,j-2) +  (16/12)*M(:,j-1) -   (2.5 )*M(:,j)  +  (16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      end
      vdah(:,end-1) =  vdah(:,end-1) + (-1/12)*M(:,end-3)+ (16/12)*M(:,end-2) - (2.5 )*M(:,end-1)  + (16/12)*M(:,end);
      vdah(:,end) =     vdah(:,end) +  (-1/12)*M(:,end-2) + (16/12)*M(:,end-1) - (2.5  - 1/12 )*M(:,end);
        