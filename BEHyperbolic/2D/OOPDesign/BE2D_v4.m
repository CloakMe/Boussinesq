function [tt, max_v, t, vu]=BE2D_v4(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0)

delete F:\MATLAB6p5\work\MT\2D\SOL\*

beta = beta1/beta2;
sy = length(y);
sx = length(x);

max_v(1) = max(max(u_t0));
   
  % dhb=deltah2d(sx,sy,1);
    W = eigvect(sx);
    DD = zeros(1,sx);
    IDH = DD;
    %I = eye(sx);
 
    for i = 1:sx
        %DD(i) = W(:,i)'*BMM([1 -4 1],W(:,i)')';%dhb*W(:,i);
        IDH(i) = W(:,i)'*BMM([(-1) (h^2+4)  (-1)],W(:,i)')';%dhb*W(:,i);
    end
         
    vdah = zeros(sx,sy);
    b = vdah;
    VV = vdah;
    vz = vdah; vmo = vdah; vu = vdah;
    

    [d2vz, d3vz, d4vz] = calc_der2d(u_t0,dudt_t0,W,IDH,3,vdah,h,sx,al,beta,4);
 
    v2 = u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz;
    
    %figure(2)
    %mesh(x,y,(tau^2/2)*d2vz')
    %title('(tau^2/2)*d^2v/dt^2 , t=tau')
    %xlabel('x');            ylabel('y');
    clear('dudt_t0'); clear('d2vz');
    %figure(2)
    %mesh(x,y,(dh_gmI(u_t0,0,vdah,vdah))/h^2');
    %title('DH * v')
    %xlabel('x');            ylabel('y');
    
    e=1;
    t(1)=0;t(2)=tau;
    k=2;
    vmo = u_t0;
    vz = v2;
    clear('u_t0'); clear('v2');
    vdah = zeros(sx,sy);
    while(t(k)<t_end)

        %vzl = [vz(2:end); 0];
        %vzr = [0; vz(1:end-1)];
        %g1 = (bt*al*(vzl.^2 + 10*vz.^2 + vzr.^2) + (bt-1)*(vzl + 10*vz +  vzr))/12;
        vu=beta*al*vz.*vz + (beta-1)*vz;
        b = W'*(dh_gmI(vu,0,vdah)); %num_ap(vu,vdah)
    
        for j=1:sx
           VV(j,:) = diag3solv([(-1) IDH(j) (-1)],b(j,:));
        end
         
        vu = 2*vz - vmo + (tau^2/h^2)*(dh_gmI(vz,0,vdah)) + tau^2*W*VV;
        
        IS_NaN = max(max(isnan(vu)));
        max_v(k) = max(max(vu));
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            %figure(2)
            %mesh(x,y,(dh_gmI(vz,0,vdah))');
            %title('DH * v')
            %xlabel('x');            ylabel('y');
            %figure(4)
            %mesh(x,y,((vu - vmo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            
            str = ['SOL\MM_' num2str(tt(e)) '.mat'];
            save (str,  'vu');%(1:xstep:end,1:xstep:end);
            %II(e)=sum(vz)*h;
            %EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            %E(e) = EL(e) + NLE(vz,vu,sx,alpha,beta);
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
            vmo = vz;
            vz  = vu;
    end
  
    
   clear('W');    clear('vmo'); 
    sol_size = size(vu)

  
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
        