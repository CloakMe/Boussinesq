clear;clc;
% constants
x_st =  -42.0;   y_st = -42.0;
x_end = 42.0;   y_end = 42.0;
pw = 1;%625
tau = 0.1; h = 0.2; x=x_st:h:x_end; y=y_st:h:y_end; t_end=8.4;
al = -1;%99979 izb
beta1 = 1;beta2 = 1;
c = 0.3;
c1 = 0.0;
vc=1;
ord=4; %O( h^ord + tau^ord )

estep = max(floor((1/tau)/10),1); %zapazwat se 20 stypki za edinitsa vreme
xstep = max(floor((1/h)/1),1);
sy = length(y)
sx = length(x)

  % HA4AJLHu YcJLOBuQ _____________
   ox1 = ones(1,sy);
    X = x'*ox1;
   if(sx~=sy)
        oy1 = ones(1,sx);
        Y = (y'*oy1)';
    else
        Y = (y'*ox1)';
   end
   
   %load 'sol_v.mat' v;
   %IC = [-2.39408114000489091743872904 0]; cut = 13.2;%   al = 0.99751;   beta = 3; alpha = 1;  cut = 13.2;
   %IC = [1.1970334615777699 0];  cut = 50;     %al = -1.994;   beta = 1/2; alpha = -2.0; 
   %u_t0 = pola2cart(x,y,IC,cut);dudt_t0 = 0*u_t0;
   if( vc == 2 )
       %IC = [2.3943685931871991900975 0];  c=0;     al = -1.0;   beta = 3; cut =10.4; gg35=@be100;   
       %IC = [2.19456515743190183620470175 0];  cut = 9.50; %  al = -1.0;   beta = 3; c=0.5
       %[u_t0 ] = pola2cart_v3(x,y,c,IC,al,beta,cut,gg35);
       
       IC = [2.29842731980465031504 0];  c=0.0;     al = -1.0;   beta = 3; cut =9.5; gg35=@ch1000;
       [u_t0] = pola2cart_v5(x, y, beta,al*beta, X, Y);
       vdah2 = zeros(size(u_t0));
       dudt_t0 = 0*u_t0;
       %dudt_t0 = (-sqrt(beta)*c/h) * (  dMdy(u_t0,vdah2)  );
   else
       if(vc==1)
           u_t0 = u_ex2d_mat_vc(X,Y,c,beta1);
           if(c==0)
            dudt_t0 = 0*u_t0;
           else
               tic 
               dudt_t0 =dudt2d_mat3_vc_v2(X,Y,c,beta1);
               toc
           end
       else
          if(vc==0)
            u_t0 = u_ex2d_mat_v2(X,Y,c,beta1);
            if(c==0)
               dudt_t0 = 0*u_t0;
            else
               tic 
               dudt_t0 =dudt2d_mat3_v2(X,Y,c,beta1);
               toc
            end
          end
       end
   end
  
   clear('X');    clear('Y');
   %figure(2)
   %mesh(x,y,u_t0');
   %mesh(x,y,dudt_t0');
   %title('IC')
   %xlabel('x');            ylabel('y');
   
   if(mod(x(end),h) == 0 && mod(y(end),h) == 0)
       [yo,indx] =min(abs(x));
       [yo,indy] =min(abs(y));
       dudt_t0(indx,indy);
       dudt_t0(indx,indy) = 0.000;
   end
   IS_dudt_NaN = max(max(isnan(dudt_t0)))
   IS_u_NaN = max(max(isnan(u_t0)))

   dscrtParams = BEDiscretizationParameters( x, y ,h, ord ,tau ,t_end, estep );
   eqParams = BEEquationParameters( al, beta1, beta2, c );
   ic = BEInitialCondition( u_t0 ,dudt_t0, 0.01 );
   
   engine = BEEngine( dscrtParams, eqParams, ic );
   % _____________________________________
   tic
  %(VS) vector scheme: -->  O(tau + h^2)
  %[tt, max_vv, t ,v1l, v2l]  = BE2D_v6(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0); ver = 6;
  %(VC) Explicit method with variable change applied -->  O(tau^2 + h^2)  tau<function(h,beta)<h ..
  %[tt, max_v, t, vl]  = BE2D_v4(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 4;
  %(NVC) Explicit method NO variable change -->  O(tau^2 + h^2) 
  %[tt, max_v, t, vl]  = BE2D_v3(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 3;
  %Taylor method variable change applied --> O(tau^4 + h^2)  tau<function(h,beta)<h ..
  %[tt, max_v, t, vl]  = BE2D_t1(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 1;
  % Taylor method variable change applied --> O(tau^4 + h^4)  tau<function(h,beta)<h ..
  [tt, max_v, t, EN, II, vl, dvl] = engine.BE2D_t2( );  ver = 2;
  % Taylor method variable change applied --> O(tau^4/tau^4 + h^8)  tau<function(h,beta)<h ..
  %[tt, max_v, t, EN, II, vl, dvl]  =  BE2D_t8(x,y,h,tau,t_end,beta1,beta2,al,c,c1,ord,0,estep,u_t0,dudt_t0);  ver = 2;
  % Taylor method variable change applied --> O(tau^4 + h^4)  tau<function(h,beta)<h ..
  %[tt, max_v, t, vl, dvl]  =  CH2D_t2(x,y,h,tau,t_end,beta1,beta2,al,c^2,ord,estep,u_t0,dudt_t0);  ver = 222;
  %(NVC sit) Explicit method NO variable change -->  O(tau^2 + h^2) 
  %[tt, max_v, t, vl]  = BE2D_v3_sit(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 33;
  %(VC sit) Explicit method with variable change applied -->  O(tau^2 + h^2)  tau<function(h,beta)<h ..
  %[tt, max_v, t, EN, II, vl]  = BE2D_v4_sit(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 44;
 
 toc
   for k = 1:sx
       if( -10^(-11)<x(k) )
           ij=k;break;
       end
   end
   for k = 1:sy
       if( -10^(-11)<y(k))
           lo=k;break;
       end
   end
 if(length(x)>floor(10/h))
    mag = floor(10/h);  %mag2 = floor(10/h);
    xx=x(ij-mag+1:ij+mag-1); yy=y(lo-mag+1:lo+mag-1);  
    if( ver == 6)
        vv1l = v1l(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag-1); 
        vv2l = v2l(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag-1);
    else
        vvl = vl(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag-1); 
    end
 else
     mag = 0;xx=x;yy=y;
 end
    maxv =  max(max(abs(u_t0)));
    if( ver == 6)
        figure(2)
        mesh(xx,yy,vv1l')
        title('solution');
        xlabel('x');            ylabel('y');
        figure(3)
        mesh(xx,yy,vv2l')
        title('solution');
        xlabel('x');            ylabel('y');
        figure(4)
        hold on;
        plot(t(1:end-1),max_vv(1,:),'b',t(1:end-1),max_vv(2,:),'r')
        hold off;
        title('Evolution of the maximum ');
        xlabel('time "t"');  ylabel('max(v)');
        return;
    else
        figure(3)
        mesh(xx,yy,vvl')
        title('solution');
        xlabel('x');            ylabel('y');
        figure(4)
        %hold on;
        plot(t(1:end-1),max_v,'k')
        %hold off;
        title('Evolution of the maximum');
        xlabel('time "t"');  ylabel('max(v)');
        figure(5)
        %hold on;
        plot(tt,EN,'k')
        %hold off;
        title('Energy functional');
        xlabel('time "t"');  ylabel('EN');
    end
    
    return;
    mag2=mag;
    figure(6)
    clear('F')
    camva('manual');
         for j = 1:size(tt,2)
            if(ver == 6)
                 load(['SOL\MM2_' num2str(tt(j)) '.mat']);
                 %vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
                 vv2u = v2u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1);
                 mesh(xx,yy,vv2u');
                 %mesh(x(1:xstep:end),y(1:xstep:end),v2u(1:xstep:end,1:xstep:end)');
                 
            else             
                load(['SOL\MM_' num2str(tt(j)) '.mat']);
                vvu = vu(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
                mesh(xx,yy,vvu');
                %mesh(x(1:xstep:end),y(1:xstep:end),vu(1:xstep:end,1:xstep:end)');
                
            end
            %mesh(x,y,vu');
            title(['Solution at time: ',num2str(tt(j))]);
            xlabel('x'); ylabel('y');
            view(87,50);
            %view(1,90);
            axis([xx(1) xx(end) yy(1) yy(end) -.7, 2.5]);
            colorbar;
            caxis([-.7, 1.7]);
            F(:,j) = getframe;
            j=j+1;
            clear('v2u');clear('vu');
         end
     movie(F,1,15)
     return;
     
     
%avi movie:
fig_ss=figure(6);
clear('F')
%set(gcf, 'OuterPosition', [0.0 30.0 1200.0 950.0]);
set(fig_ss,'DoubleBuffer','on');
%set(gca,'replace','Visible','off');
winsize = get(fig_ss,'Position');
mov = avifile('example5.avi');
mov.Quality = 100;
mov.FPS = 5;
    camva('manual');
         for j = 1:size(tt,2)
          if(tt(j)>20) break; end
            if(ver == 6)
                 load(['SOL\MM2_' num2str(tt(j)) '.mat']);
                 vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag-1); 
                 vv2u = v2u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag-1);
                 mesh(xx,yy,vv2u');
                 %mesh(x(1:xstep:end),y(1:xstep:end),v2u(1:xstep:end,1:xstep:end)');
            else             
                load(['SOL\MM_' num2str(tt(j)) '.mat']);
                vvu = vu(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
                mesh(xx,yy,vvu');
                %mesh(x(1:xstep:end),y(1:xstep:end),vu(1:xstep:end,1:xstep:end)');
                
            end
            maxv=max(max(vvu)); minv=min(min(vvu));
            %mesh(x,y,vu');
            title(['Solution at time: ',num2str(tt(j))]);
            xlabel('x'); ylabel('y');
            view(90,90);
            %view(1,90);
            axis([xx(1) xx(end) yy(1) yy(end) -.5, 2.36]);
            colorbar;
            caxis([-.5, 2.36]);
            %F(:,j) = getframe;
            G = getframe(fig_ss,[18,4,508 420]);
            j=j+1;
            clear('v2u');clear('vu');
            mov = addframe(mov,G);
         end
         mov = close(mov);
     
     movie(mov,1,20)
%     return;
     load(['SOL\MM_' num2str(size(tt,2)) '.mat']);
     figure(7)
         for j = 1:size(vu,2)
            %['SOL\MM_' num2str(j) '_mat']
            
            plot(x,vu(:,j));
            xlabel('x'); ylabel('z');
            axis([x_st x_end -3.1 maxv]);
            F(:,j) = getframe;
            j=j+1;
         end
         clear('vu');
     movie(F,1,1)
     
     load(['SOL\MM_' num2str(size(tt,2)) '.mat']);
     figure(7)
         for j = 1:size(vu,1)
             %['SOL\MM_' num2str(j) '_mat']
            
            plot(y,vu(j,:));
            xlabel('y'); ylabel('z');
            axis([x_st x_end -3.1 maxv]);
            F(:,j) = getframe;
            j=j+1;
         end
         clear('vu');
     movie(F,1,1)
 
  figure(7)
  clear('F')
         for j = 1:size(tt,2)
             %['SOL\MM_' num2str(j) '_mat']
            load(['SOL\MM_' num2str(j) '.mat']);
            %plot(x,vu(:,(1+end)/2)');
            %xlabel('x'); ylabel('v');
            plot(y,vu((1+end)/2,:)');
            xlabel('y'); ylabel('v');
            axis([x_st x_end  -0.5 maxv]);
            F(:,j) = getframe;
            j=j+1;
            clear('vu');
         end
     movie(F,1,5)
     
     
     figure(7)
         %mesh(x(268:end-267),y(268:end-267),u_t0(268:end-267,268:end-267)');
         plot(y(268:end-267),u_t0((1+end)/2,268:end-267)','.k');
    
        for j = 2:2:12
             %['SOL\MM_' num2str(j) '_mat']
             figure(j/2)
            load(['SOL_VS_3\MM_' num2str(j) '.mat']);
            load(['SOL_VS_3\MM2_' num2str(j) '.mat']);
            %load(['SOL_NVC_3\MM_' num2str(j) '.mat']);
                  mean_vs=(v1u+v2u)/2;
            mesh(x,y,vu'-mean_vs');
            MAXd = max(max(v2u'));%-mean_vs'));
            MINd = min(min(v2u'));%-mean_vs')); 
            %mesh(x,y,vu');
            %title(['Difference at time: ',num2str(j)]);
            xlabel('x'); ylabel('y');
            view(-90,90);
            %view(1,90);
            axis([x_st x_end y_st y_end -.5 2.4]);
            colorbar;
            caxis([-.5 2.4]);
            
            %clear('vu');clear('v1u');clear('v2u');
        end
        
         figure(7)
         %mesh(x(268:end-267),y(268:end-267),u_t0(268:end-267,268:end-267)');
         plot(y(268:end-267),u_t0((1+end)/2,268:end-267)','.k');
    hold on;
    yoyo = ['-k' '.k' '-k' '.k' '-k' '.k']; 
        for j = 2:2:12
            %['SOL\MM_' num2str(j) '_mat']
             
            %load(['SOL_VS\MM_' num2str(j) '.mat']);
            load(['SOL_VS\MM2_' num2str(j) '.mat']);
            %load(['SOL_NVC\MM_' num2str(j) '.mat']);
           % figure(j/2)
            %mean_vs=(v1u+v2u)/2;
            %mesh(x(266:end-266),y(266:end-266),v2u(266:end-266,266:end-266
            %)'); %vu - v2u
            plot(y(268:end-267),v2u((1+end)/2 ,268:end-267)',yoyo(j-1:j))
            MAXd = max(max(v2u'));
            MINd = min(min(v2u')); 
            %mesh(x,y,vu');
            title(['Cross-section at x=0 ']);
            %xlabel('x'); ylabel('y');
            xlabel('y'); ylabel('v');
            %view(124,24);
            %view(124,24);
            %axis([x(266) x(end-266) y(266) y(end-266) -0.2 2.4]);
            %colorbar;
            %caxis([-0.2, 2.4]);
            
            %clear('vu');clear('v1u');clear('v2u');
        end
         %   [tt, max_vv, t ,v1l, v2l] 
         tt2=tt;
         max_vv2=max_vv;
         t2=t;
         v1l2=v1l;
         qs = (sx-1)/4;
         v2u = v2l;
        
        size(v2l(1+qs:end-qs,1+qs:end-qs)')
        size(v2u')
        figure(2)
        mesh(x(1+qs+66:end-qs-66),y(1+qs+66:end-qs-66),v2l(1+qs+66:end-qs-66,1+qs+66:end-qs-66)' - v2u(67:end-66,67:end-66)')
        title('diff');
        xlabel('x');            ylabel('y');
        figure(3)
        
        figure(4)
        plot(t(1:end-1),max_vv(1,:),'--k',t(1:end-1),max_vv(2,:),'k')
        title('Evolution of the maximum ');
        xlabel('time "t"');  ylabel('max(v)');