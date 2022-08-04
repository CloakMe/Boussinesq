function [tt, max_v, t, vu]=BE2D_v3(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0)

delete F:\MATLAB6p5\work\MT\2D\SOL\*

 %x=x_st:h:x_end; y=y_st:h:y_end; 
  
sy = length(y);
sx = length(x);

%sx2 = sx*sy;
max_v(1) = max(max(u_t0));
      
    W = eigvect(sx);
    DD = zeros(1,sx);
    IDH = DD;
    %I = eye(sx);
 
    for i = 1:sx
        %DD(i) = W(:,i)'*BMM([1 -4 1],W(:,i)')';%dhb*W(:,i);
        IDH(i) = W(:,i)'*BMM([((-1)) (h^2/beta1+4)  ((-1))],W(:,i)')';%dhb*W(:,i);
    end
         
    vdah = zeros(sx,sy);
    b = vdah;
    VV = vdah;
    vz = vdah; vmo = vdah; vu = vdah;
        
    %{
    vz=(u_t0 + al*u_t0.*u_t0 + (beta2/h^2)*(-dh_gmI(u_t0,1/beta2,vdah,vdah)))/beta1;% + (beta1-1)*u_t0;
    
    vdah = zeros(sx,sy);
    vdah(:,1:end-1) = vz(:,2:end);    vdah(:,2:end) = vdah(:,2:end) + vz(:,1:end-1);
   
   % tic;b = diag(DD)*W'*vz + W'*vdah;toc;
    
    
    tic;b = W'*(dhb*vz + vdah);toc;
            
    for j=1:sx
        VV(j,:) = diag3solv([(-1) IDH(j) (-1)],b(j,:));
    end

    dv2 = W*VV;
    v2 = u_t0 + tau*dudt_t0 + (tau^2/2)*dv2;
    
    %}
    
    vdah = zeros(sx,sy);
 
    [d2vz, d3vz, d4vz] = calc_der2d_v2(u_t0,dudt_t0,W,IDH,3,vdah,h,sx,al,beta1,beta2,4);
  
    v2 = u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz;
   %{
   figure(2)
    mesh(x,y,(tau^2/2)*d2vz')
    title('(tau^2/2)*d^2v/dt^2 , t=tau')
    xlabel('x');            ylabel('y');
     %}   
    e=1;
    t(1)=0;tt(1)=tau;t(2)=tau;
    k=2;
    vmo = u_t0;
    vz = v2;
    clear('u_t0'); clear('v2');
    clear('dudt_t0'); clear('d2vz');
 
    while(t(k)<t_end)
       
        vu = al*vz.*vz;
        vu = vu + (beta2/h^2)*(-dh_gmI(vz,h^2/beta2,vdah));% + (beta1-1)*u_t0;
           
        b = W'*(dh_gmI(vu,0,vdah))/beta1;
            
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],b(j,:));
        end
    
        vu = 2*vz - vmo + tau^2*W*VV;
        
        IS_dudt_NaN = max(max(isnan(vu)));
        max_v(k) = max(max(vu));
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            %{
            figure(2)
            mesh(x,y,(dh_gmI(vu,0,vdah)/beta1)');
            title('DH * v')
            xlabel('x');            ylabel('y');
            figure(4)
            mesh(x,y,((vu - vmo)/(2*tau))');
            title('v_t')
            xlabel('x');            ylabel('y');
            colorbar;
            %}
            str = ['SOL\MM_' num2str(tt(e)) '.mat'];
            save (str,  'vu');%(1:xstep:end,1:xstep:end);
            %II(e)=sum(vz)*h;
            %EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            %E(e) = EL(e) + NLE(vz,vu,sx,alpha,beta);
            e=e+1;
        end
        if(IS_dudt_NaN == 1)
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
  
