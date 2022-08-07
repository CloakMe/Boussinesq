function [tt, max_v, t, v1u, v2u]=BE2D_v6(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0)


delete SOL\*
  
sy = length(y);
sx = length(x);

    d2x_s = deltah2d_xy(7,h);
    d2x_ss = d2x_s*d2x_s;
    eye7 = eye(7);
    sqD = sqrt(  ( 1 + beta1/tau^2)^2 - 4*beta2/tau^2 );
    X1 =  (1 + beta1/tau^2 + sqD)/( 2 * h^2 / tau^2);
    X2 =  (1 + beta1/tau^2 - sqD)/( 2 * h^2 / tau^2);
    %(beta2/h^2)*(d2x_s-X1*eye7)*(d2x_s-X2*eye7)
    A1 = (-d2x_s+X1*eye7);
    A2 = (-d2x_s+X2*eye7);
    A = h^2*eye7/tau^2 - d2x_s + beta2*d2x_ss/h^2 - beta1*d2x_s/tau^2;
    sA = A(3,1:5);
    sA1 = A1(2,1:3);
    sA2 = A2(2,1:3);
    coeff = (h^2/beta2);
    a11 = A(1,1);
    
    %Iy = eye(sy);
    W = eigvect(sx);
    IDH = zeros(1,sx);
    for i = 1:sx
       IDH(i) = W(:,i)'*BMM([(-1) (h^2/beta1+4)  (-1)],W(:,i)')';
    end
    
    max_v(1,1) = max(max(abs(u_t0)));
    max_v(2,1) = max(max(abs(u_t0)));
   
    %[D2x]=deltah2d_xy(sx,h);
    %[D2y]=deltah2d_xy(sy,h);
         
    VV = zeros(sx,sy);
    vdah = zeros(sx,sy);
    %tic
    %D2x*u_t0;
    %u_t0*D2y;       
    %toc
    %tic
    %diff2x(u_t0,vdah);
    %diff2y(u_t0,vdah);
    %toc
    b = VV;      b2 = VV;
    v1z = VV; v1mo = VV; v1u = VV;
    v2z = VV; v2mo = VV; v2mt = VV; v2u = VV;        
  
    [d2vz, d3vz, d4vz] = calc_der2d_v2(u_t0,dudt_t0,W,IDH,3,vdah,h,sx,al,beta1,beta2,4);

    v1z = u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz;
    
    v2z = v1z;
    v1mo = u_t0; v2mo = u_t0;
    v2mt =  u_t0 - tau*dudt_t0 + (tau^2/2)*d2vz - (tau^3/6)*d3vz + (tau^4/24)*d4vz;
    v1mt = v2mt;
    %figure(1)
    %mesh(x,y,v1z')
    clear('W'); clear('u_t0'); clear('dudt_t0');
    
    e=1;
    t(1)=0;tt(1)=tau;t(2)=tau;
    k=2;
    %{
    clear('rt');tav = .08;rt(1)=0;
    op=1;
   for o = 2:201
       op=op+1;
      rt(op) = (op-1)*tav;
   end
   %}
    while(t(k)<t_end-1.0e-08)

        b2 = (2*v1z - v1mo)/tau^2;
        VV = diff2x ( al* v1z.^2 , vdah );
        d2vz = diff2y( al* v2z.^2 , vdah  );
        
        gm = (h^2/beta2)*( 1 + beta1/tau^2);
        
        d3vz = diff2x ( -(beta2/h^2) * diff2y(v2z+v1z , vdah) - beta1*b2 , vdah );
        d4vz = diff2y( (-beta2/h^2) * diff2y(v2z , vdah) + v2z + (beta1/tau^2)*(v2z - 2*v2mo + v2mt) , vdah  );
        
        b = b2*h^2 + VV + d2vz + d3vz + d4vz;

       f1=0; f2=0;
        for j=1:(sy+1)/2
            if( norm( b(:,j) - b(:,end - j + 1), 2) < 1.00e-10 )
                %v1u(:,j) = diag3solv(sA1,b(:,j));  v1u(:,j) = coeff*diag3solv(sA2,v1u(:,j));
                v1u(:,j) = pentsolv(a11,sA,b(:,j));
                v1u(:,end - j + 1) = v1u(:,j);%pentsolv(a11,sA,b(:,end - j + 1));
                f1 = f1+1;
            else
                v1u(:,j) = pentsolv(a11,sA,b(:,j));
                v1u(:,end - j + 1) = pentsolv(a11,sA,b(:,end - j + 1));
                
                %v1u(:,j) = diag3solv(sA1,b(:,j));  v1u(:,j) = coeff*diag3solv(sA2,v1u(:,j));
                %v1u(:,end - j + 1) = diag3solv(sA1,b(:,end - j + 1));
                %v1u(:,end - j + 1) = coeff*diag3solv(sA2,v1u(:,end - j + 1));
                f2 = f2+1;
            end
        end

        f1=0; f2=0;
        gm = (h^2/beta2)*( 1 + beta1/tau^2);
        b = (2*v2z - v2mo)/tau^2;
        
        %d2vz = al*diff2y((v2z.^2)  , vdah );
        %VV = diff2x (  al* v1z.^2 , vdah );
        
        d3vz = diff2x( ( -beta2/h^2) * diff2x(v1u , vdah) + (1 + beta1/tau^2)*v1u -( beta2/h^2)* diff2y(v2z+v1u , vdah)   - beta1*b2   , vdah );
        d4vz = diff2y(  -beta1 * b , vdah );
        b2 = b*h^2  + VV +  d2vz + d3vz + d4vz;
  
        for i=1:(sx+1)/2
            if( norm(b2(i,:) - b2(end - i + 1,:)  , 2) < 1.00e-10 )
                v2u(i,:) = pentsolv(a11,sA,b2(i,:));
                %v2u(i,:) = diag3solv(sA1,b2(i,:));    v2u(i,:) = coeff*diag3solv(sA2,v2u(i,:));
                v2u(end - i + 1,:) = v2u(i,:);%pentsolv(a11,sA,b2(end - i + 1,:));
                f1 = f1+1;
            else
                v2u(i,:) = pentsolv(a11,sA,b2(i,:));
                v2u(end - i + 1,:) = pentsolv(a11,sA,b2(end - i + 1,:));
                
                %v2u(i,:) = diag3solv(sA1,b2(i,:));    v2u(i,:) = coeff*diag3solv(sA2,v2u(i,:));
                %v2u(end - i + 1,:) = diag3solv(sA1,b2(end - i + 1,:));
                %v2u(end - i + 1,:) = coeff*diag3solv(sA2,v2u(end - i + 1,:));
                f2 = f2+1;
            end
        end
                     
        IS_NaN = max(max(isnan(v1u))) + max(max(isnan(v2u)));
        max_v(1,k) = max(max(abs(v1u)));
        max_v(2,k) = max(max(abs(v2u)));
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            %figure(2)
            %mesh(x,y,(vdah + dhb*vz)');
            %title('DH * v')
            %xlabel('x');            ylabel('y');
            %figure(5)
            %mesh(x,y,v1u');
            %mesh(x,y,((v1u - v1mo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            %figure(6)
            %mesh(x,y,v2u');
            %mesh(x,y,((v2u - v2mo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            
            str = ['SOL\MM_' num2str(tt(e)) '.mat'];
            save (str,  'v1u');%(1:xstep:end,1:xstep:end);
            str = ['SOL\MM2_' num2str(tt(e)) '.mat'];
            save (str,  'v2u');%(1:xstep:end,1:xstep:end);
            %II(e)=sum(vz)*h;
            %EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            %E(e) = EL(e) + NLE(vz,vu,sx,alpha,beta);
            %if(abs(tt(e)-t_end) < 1.00e-010)
            %       tt(end)
            %       t(end)
            %end
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
        if(IS_NaN > 0)
           warning('ERROR; BE2D NaN values; Stopping! ');
           t(end) 
           return;
        end
          k=k+1;  
           
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k)))==0)
                yoyo=t(k)
            end
            v1mo = v1z;
            v1z  = v1u;
            
            v2mt = v2mo;
            v1mt = v1mo;
            v2mo = v2z;
            v2z  = v2u;
    end
    
   clear('v1mo');  clear('v2mo'); 
    sol_size = size(v1u)

function vdah=diff2y(M,vdah)
      j=1;
      vdah(:,j) =  - (2)*M(:,j) +  M(:,j+1);
      for j=2:size(M,1)-1
          vdah(:,j) =  +  M(:,j-1) - (2)*M(:,j) +M(:,j+1);
      end
      j=size(M,1);
      vdah(:,j) =  +M(:,j-1) - (2)*M(:,j) ;
      %vdah = vdah/h^2;
      
function vdah=diff2x(M,vdah)
      j=1;
      vdah(j,:) = - (2)*M(j,:) +  M(j+1,:);
      for j=2:size(M,1)-1
          vdah(j,:) =  M(j-1,:) - (2)*M(j,:) +  M(j+1,:);
      end
      j=size(M,1);
      vdah(j,:) =  M(j-1,:) - (2)*M(j,:) ;
      %vdah = vdah/h^2;
      
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
      
function resmm=d2y_gmI(M,gm,vdah)
      vdah(:,1:end-1) = M(:,2:end);    
      vdah(:,2:end) = vdah(:,2:end) + M(:,1:end-1);
      resmm = vdah - (2 +gm)*M;
      
function resmm=d2x_gmI(M,gm,vdah)
      vdah(1:end-1,:) = M(2:end,:); 
      vdah(2:end,:) = vdah(2:end,:) + M(1:end-1,:);
      resmm = vdah- (2 +gm)*M;
        