clear;clc;
% constants
% st ; end
%-20 ; 15
start_x=-100; end_x = 100;

for j = 1:4
clear('x');clear('uex2');clear('u_t0');clear('dudt_t0');clear('v');clear('dtv');clear('tt');
clear('dh');clear('sdh');clear('Idh');clear('sIdh');
h = 0.2/2^(j-1);  tau = 0.2/2^(j-1);  x = start_x:h:end_x;  t_end=10;
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
  %[v,tt,E,II] = BE1D_vnum1(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
  
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
    [v,dtv,tt,E,II] = BE1D_tv2(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
    
    dh = deltaOh4(7);
    sdh = dh(3,1:5);
    sdh11 = dh(1,1);
    Idh = h^2*eye(7) - deltaOh4(7);
    sIdh = Idh(3,1:5);
    sIdh11 = Idh(1,1);
    
    [d2v1, d3v1, d4v1, d5v1] = calc_der(v,dtv,sdh,sIdh,sdh11,sIdh11,h,alpha,beta);
       
     VL2.(genvarname(['v' num2str(j)])) = d2v1;
     VL3.(genvarname(['v' num2str(j)])) = d3v1;
     VL4.(genvarname(['v' num2str(j)])) = d4v1;
     VL5.(genvarname(['v' num2str(j)])) = d5v1;
   %Integral = II(end)
    %Energy = E(end)
    uex2 =  u_ex(x+shift,tt(end),c,alpha,beta1,beta2);% + u_ex(x(l)-5,t_end,-1.5);
   
    Error = max(abs(uex2' - v(:,end)))
end

if(length(VL2.v1)~= length(VL2.v2))
   if(2*length(VL2.v1)-1== length(VL2.v2))
    res11 = uex2(1:8:end)' - VL2.v1(:);
    res12 = uex2(1:8:end)' - VL2.v2(1:2:end);
    res21 = uex2(1:4:end)' - VL2.v2(:);
    res22 = uex2(1:4:end)' - VL2.v3(1:2:end);
 
    E1_1 = norm(res11(:),2);
    E1_2 = norm(res12(:),2);
    E2_1 = norm(res21(:),2);
    E2_2 = norm(res22(:),2);

    conv_1 = (log(E1_1)-log(E1_2))/log(2)
    conv_2 = (log(E2_1)-log(E2_2))/log(2)
    if(j==4)
        res31 = uex2(1:2:end)' - VL2.v3(:);
        res32 = uex2(1:2:end)' - VL2.v4(1:2:end);
        
        E3_1 = norm(res31(:),2);
        E3_2 = norm(res32(:),2);
        
        conv_3 = (log(E3_1)-log(E3_2))/log(2)
    end
   end
else

    res1 = uex2' - VL2.v1(:);
    figure(1)
    plot(x,res1)
    res2 = uex2' - VL2.v2(:);
    figure(2)
    plot(x,res2)
    res3 = uex2' - VL2.v3(:);
    figure(3)
    plot(x,res3)
    E1_1 = norm(res1(:),2);
    E2_1 = norm(res2(:),2);
    E3_1 = norm(res3(:),2);

    conv_1 = (log(E1_1)-log(E2_1))/log(2)
    conv_2 = (log(E2_1)-log(E3_1))/log(2)
    if(j==4)
        res4 = uex2' - VL2.v4(:);
        E4_1 = norm(res4(:),2);
        
        conv_3 = (log(E3_1)-log(E4_1))/log(2)
    end
end

%{
    res1 = VL5.v1(:) - VL5.v2(1:2:end);
    res2 = VL5.v2(1:2:end) - VL5.v3(1:4:end);
    E1_1 = norm(res1(:),2)
    E2_1 = norm(res2(:),2)

    conv = (log(E1_1)-log(E2_1))/log(2)
        if(j==4)
            re1 = VL5.v2(:) - VL5.v3(1:2:end);
            re2 = VL5.v3(1:2:end) - VL5.v4(1:4:end);
            E1_2 = norm(re1(:),2)
            E2_2 = norm(re2(:),2)

            conv = (log(E1_2)-log(E2_2))/log(2)
        end

%}
    