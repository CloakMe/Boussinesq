clear;
L=40;  
T=5;  

h=0.2; tau=0.1;
[x1,t1,u1] = BPE_symplectic_periodic_condition( L, T, h, tau);

h=0.1; tau=0.05;
[x2,t2,u2] = BPE_symplectic_periodic_condition( L, T, h, tau);

h=0.05; tau=0.025;
[x3,t3,u3] = BPE_symplectic_periodic_condition( L, T, h, tau);

res_1 = u1 - u2(1:2:end);
res_2 = u2 - u3(1:2:end);

norm0402_L2 = 4*h*norm(res_1(:),2);
norm0201_L2 = 2*h*norm(res_2(:),2);
fprintf('Solution Convergence:\n');

fprintf('||v_04 - v_02||_L2 = %.6f \n', norm0402_L2);
fprintf('||v_02 - v_01||_L2 = %.6f \n', norm0201_L2); %%.8e
conv_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2);
fprintf('Conv_L2 = %.8e \n\n', conv_L2);

norm0402_Inf = max(max(abs(res_1(:))));
norm0201_Inf = max(max(abs(res_2(:))));
fprintf('||v_04 - v_02||_Inf = %.6f \n', norm0402_Inf);
fprintf('||v_02 - v_01||_Inf = %.6f \n', norm0201_Inf);
conv_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2);
fprintf('Conv_Inf = %.8e \n', conv_Inf);
