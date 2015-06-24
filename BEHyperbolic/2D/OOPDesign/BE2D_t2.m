function [tt, max_v, t, EN, II, vu, dtv]=BE2D_t2( dscrtParams, eqParams, ic )
                                                %( dscrtParams, eqParams, ic )
                                                %(x,y,h,tau,t_end,beta1,beta2,al,ord,estep,u_t0,dudt_t0)                                                

delete SOL\*

beta = beta1/beta2;
sy = length(y);
sx = length(x);

%sx2 = sx*sy;
max_v(1) = max(max(abs(u_t0)));
    [ff,dhb]=deltah2d_xy(sx,h);
    [W,w] = eig(-dhb);
    s_idh = 5;
    DD = zeros(1,sx);
    IDH = DD;
    %I = eye(sx);

    for i = 1:sx
        %DD(i) = W(:,i)'*BMM([1 -4 1],W(:,i)')';%dhb*W(:,i);
        IDH(i) = W(:,i)'*BMM([(1/12) (-16/12) (h^2+5)  (-16/12) (1/12)],W(:,i)',(h^2+59/12))';%dhb*W(:,i);
        IDHo(i) = W(:,i)'*BMM([(1/12) (-16/12) (5)  (-16/12) (1/12)],W(:,i)',(59/12))';
    end
    
        
    vdah = zeros(sx,sy);
    b = vdah;
    VV = vdah;
    vu = vdah;
        
    [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(u_t0,dudt_t0,W,IDH,s_idh,vdah,h,sx,al,beta,ord);
    if(ord == 6)
    vu = u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz + (tau^5/120)*d5vz + (tau^6/720)*d6vz;
    dtv = dudt_t0 + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;
    else
        vu =  u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz;
        dtv = dudt_t0 + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz ;
    end
        
    %figure(2)
    %mesh(x,y,(tau^2/2)*d2vz')
    %title('(tau^2/2)*d^2v/dt^2 , t=tau')
    %xlabel('x');            ylabel('y');
    clear('dudt_t0'); clear('d2vz');
        
    e=1; tt(e) = 0;
    t(1)=0;t(2)=tau;
    k=2;
    vz = vu; vmo = u_t0;
    clear('u_t0'); clear('v2');
    EN(1)=0;II(1)=0;
    while(t(k)<t_end)
        
        [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vu,dtv,W,IDH,s_idh,vdah,h,sx,al,beta,ord);
        if(ord == 6)
            vu = vu + tau*dtv + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz+ (tau^5/120)*d5vz + (tau^6/720)*d6vz;
            dtv = dtv + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz + (tau^4/24)*d5vz + (tau^5/120)*d6vz;    
        else
            vu =    vu + tau*dtv + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz;
            dtv = dtv + tau*d2vz + (tau^2/2)*d3vz + (tau^3/6)*d4vz;
        end

        IS_dudt_NaN = max(max(isnan(vu)));
        max_v(k) = max(max(abs(vu)));
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            %figure(2)
            %mesh(x,y,(dh_gmI(vz,0,vdah,vdah))');
            %title('DH * v')
            %xlabel('x');            ylabel('y');
            %figure(4)
            %mesh(x,y,((vu - vmo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            EN(e) = Eg2_vc_2d(vmo,vz,vu,IDHo,W,h,tau,al,beta,sx,sy,s_idh,vdah,vdah,vdah);
            str = ['SOL\MM_' num2str(tt(e)) '.mat'];
            save (str,  'vu');%(1:xstep:end,1:xstep:end);
            II(e)=sum(sum(vz))*h^2;
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
        vmo = vz; vz = vu; 
    end

    
   clear('W');    clear('vmo'); 
    sol_size = size(vu)
end
  
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
end
function resmm=num_ap(M,vdah)
      vdah(1:end-1,:) = M(2:end,:); 
      vdah(2:end,:) = vdah(2:end,:) + M(1:end-1,:);
      resmm = (vdah + 10*M)/12;
end