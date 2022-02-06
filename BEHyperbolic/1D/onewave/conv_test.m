clear;clc;
load v_t10_h1.mat;
load v_t10_h05.mat;
load v_t10_h025.mat;

v_h025=v_t10_h025_end;
v_h05=v_t10_h05_end;
v_h1=v_t10_h1_end;

size(v_h1)
size(v_h05)
size(v_h025)

E1 = norm(v_h1 - v_h05(1:2:end),2)
E2 = norm(v_h05(1:2:end) - v_h025(1:4:end),2)

conv = (log(E1)-log(E2))/log(2)

clear;clc;
load v_t10_h1_t1.mat;
load v_t10_h1_t05.mat;
load v_t10_h1_t025.mat;
load v_t10_h1_t0125.mat;

v_t025=v_t10_h1_t0125_end;
v_t025=v_t10_h1_t025_end;
v_t05=v_t10_h1_t05_end;
v_t1=v_t10_h1_t1_end;

size(v_t1)
size(v_t05)
size(v_t025)
size(v_t0125)

E1 = norm(v_t1 - v_t05(1:end),2)
E2 = norm(v_t05(1:end) - v_t025(1:end),2)

conv = (log(E1)-log(E2))/log(2)

