clear;clc;
% constants
% st ; end
%-20 ; 15
start_x=-100; end_x = 100;

for j = 1:4
clear('x');clear('uex2');clear('u_t0');clear('dudt_t0');clear('v');clear('tt');
h = 0.1/2^(j-1);  tau = 0.1/2^(j-1);  x = start_x:h:end_x;  t_end=20;
beta1=1.5;   beta2=0.5;  alpha=3; beta=beta1/beta2;
%sgm = 0;
%sgm=(1-h^2/(tau^2))/12;
%sgm=h^2/(12*tau^2);
%sgm=(-h^2/(tau^2))/12;
%sgm=(beta-h^2/(tau^2))/12;
sgm =  1/12;
c=2;   shift = 20;
estep = max(floor((1/tau)/10),1); %zapazwat se 20 stypki za edinitsa vreme

    %if beta1, beta2 and alpha not specified then beta1 = 1.5; beta2 = 0.5, alpha = 3.
    
    %Tochno reshenie u_ex(x,t,c,alpha,beta1,beta2)
    % t == 0
    u_t0 = u_ex(x+shift,0,c,alpha,beta1,beta2);% + u_ex(x-shift,0,-1.5);
    %Tochno reshenie - proizwodna - dudt_ex(x,t,c,beta1,beta2,alpha)
    % t == 0
    dudt_t0 =dudt_ex(x+shift,0,c,alpha,beta1,beta2);% + dudt_ex(x(l)-shift,0,-1.5);

%==========================================================================================   
% vz - v-zero layer
% NM; O(tau^2 + h^2); 
%  g_0 = alpha*beta*vz^2  +  (beta-1)*vz      
  %[v,tt,E,II] = BE1D_v10(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   
% vz - v-zero layer
% NM; O(tau^2 + h^2); -->numerov along x-axis; nonlinear term - v^2 only :
%   g_0 --> alpha*beta*(vz(x_(i-1))^2 + 10*vz(x_i)^2 + vz(x_(i+1))^2 )/12      
   %[v,tt,E,II] = BE1D_v11(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   
% IM -->Numerov app: g_2 = [ g_0(v_up) + 10*g_0(v_zero) + g_0(v_down) ] / 12  -  if sgm = (1 - h^2/tau^2)/12 ==>  O(tau^4 + h^2) 
%  [v,tt,E,II] = BE1D_vnum1(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
  
% IM --> g_1 = [ G(v_up) - G(v_down) ]. / (v_up - v_down)  -    O(tau^2 + h^2) 
  %[v,tt,E,II] = BE1D_v3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0); 
    
% IM --> g_3 =  2 * [ G( (1/2)(v_up+v_mid) ) - G( (1/2)(v_mid+v_down) ) ]. / (v_up - v_down) 
% if sgm = (1 - h^2/tau^2)/12 ==>  O(tau^2 + h^2)
   %[v,tt,E,II] = BE1D_v4(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
    
% Taylor  -  O(tau^4 + h^2)
    %[v,dtv,tt,E,II] = BE1D_taylor(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
    
% Taylor v2  -   O(tau^4 + h^4)
    %dh = [1 -2 1]/h^2 se zamenq s dh = [-1 16 -30 16 -1]/(12*h^2) i
    %podobrqwame reda na sxodimost w prostranstwenite koordinati
    [v,dtv,tt,II] = BE1D_tv2(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
  
     VL.(genvarname(['v' num2str(j)])) = v(:,end);
     Integral.(genvarname(['v' num2str(j)])) = II(end)
    %Integral = II(end)
    %Energy = E(end)
    uex2 =  u_ex(x+shift,tt(end),c,alpha,beta1,beta2);% + u_ex(x(l)-5,t_end,-1.5);
   
    Error = max(abs(uex2' - v(:,end)))
end

if(length(VL.v1)~= length(VL.v2))
   if(2*length(VL.v1)-1== length(VL.v2))
       if(j==3) jj=2; else jj=1; end
    res11 = uex2(1:8/jj:end)' - VL.v1(:);
    res12 = uex2(1:8/jj:end)' - VL.v2(1:2:end);
    res21 = uex2(1:4/jj:end)' - VL.v2(:);
    res22 = uex2(1:4/jj:end)' - VL.v3(1:2:end);
 
    E1_1 = norm(res11(:),2);
    E1_2 = norm(res12(:),2);
    E2_1 = norm(res21(:),2);
    E2_2 = norm(res22(:),2);

    conv_1 = (log(E1_1)-log(E1_2))/log(2)
    conv_2 = (log(E2_1)-log(E2_2))/log(2)
    if(j==4)
        res31 = uex2(1:2:end)' - VL.v3(:);
        res32 = uex2(1:2:end)' - VL.v4(1:2:end);
        
        E3_1 = norm(res31(:),2);
        E3_2 = norm(res32(:),2);
        
        conv_3 = (log(E3_1)-log(E3_2))/log(2)
    end
   end
else

    res1 = uex2' - VL.v1(:);
    figure(1)
    plot(x,res1)
    res2 = uex2' - VL.v2(:);
    figure(2)
    plot(x,res2)
    res3 = uex2' - VL.v3(:);
    figure(3)
    plot(x,res3)
    E1_1 = norm(res1(:),2);
    E2_1 = norm(res2(:),2);
    E3_1 = norm(res3(:),2);

    conv_1 = (log(E1_1)-log(E2_1))/log(2)
    conv_2 = (log(E2_1)-log(E3_1))/log(2)
    if(j==4)
        res4 = uex2' - VL.v4(:);
        E4_1 = norm(res4(:),2);
        
        conv_3 = (log(E3_1)-log(E4_1))/log(2)
    end
end

%{
% V sluchai che sme namalili 4 pyti stypkata po prostranstwoto
% i programata ne pokazwa izraza "conv" za shodimostta to togawa puskame
% koda po dolu:

res1 = VL.v1(:) - VL.v2(1:4:end);
    res2 = VL.v2(1:4:end) - VL.v3(1:16:end);
    E1_1 = norm(res1(:),2)
    E2_1 = norm(res2(:),2)

    conv = (log(E1_1)-log(E2_1))/log(2)
%}

%{ 
    res31 = uex2(1:4:end)' - VL.v1(:);
        res32 = uex2(1:end)' - VL.v2(1:end);
        
        E3_1 = norm(res31(:),2);
        E3_2 = norm(res32(:),2);
        
        conv_3 = (log(E3_1)-log(E3_2))/log(2)
        %}