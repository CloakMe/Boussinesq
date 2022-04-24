clear;clc;
load vta_t10_h2_t07.mat
load vta_t10_h1_t05.mat
load vta_t10_h05_t035.mat
load vta_t10_h025_t025.mat
%oad vta_t10_h025_t025

size(vta_t10_h2_t07e)
size(vta_t10_h1_t05e)
size(vta_t10_h05_t035e)
size(vta_t10_h025_t025e)
E1 = norm(vta_t10_h2_t07e - vta_t10_h1_t05e(1:2:end),2)
E2 = norm(vta_t10_h1_t05e(1:2:end) - vta_t10_h05_t035e(1:4:end),2)

conv = (log(E1)-log(E2))/log(2)

EE1 = norm(vta_t10_h1_t05e - vta_t10_h05_t035e(1:2:end),2)
EE2 = norm(vta_t10_h05_t035e(1:2:end) - vta_t10_h025_t025e(1:4:end),2)

conv2 = (log(EE1)-log(EE2))/log(2)

plot(x,vta_t10_h025_t025e,'k',x(1:2:end),vta_t10_h05_t035e,'b',x(1:4:end),vta_t10_h1_t05e,'g')
