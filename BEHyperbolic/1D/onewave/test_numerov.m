clear;clc; %wij faila do dolu!
load vn_t10_h1.mat
load vn_t10_h05.mat
load vn_t10_h025.mat
load vn_t10_h0125.mat


size(vn_t10_h1e)
size(vn_t10_h05e)
size(vn_t10_h025e)
size(vn_t10_h0125e)
E1 = norm(vn_t10_h1e - vn_t10_h05e(1:2:end),2)
E2 = norm(vn_t10_h05e(1:2:end) - vn_t10_h025e(1:4:end),2)

conv = (log(E1)-log(E2))/log(2)

EE1 = norm(vn_t10_h05e - vn_t10_h025e(1:2:end),2)
EE2 = norm(vn_t10_h025e(1:2:end) - vn_t10_h0125e(1:4:end),2)

conv2 = (log(EE1)-log(EE2))/log(2)



clear;clc;
load vta_t10_ht05.mat
load vta_t10_ht025.mat
load vta_t10_ht0125.mat
load vta_t10_ht06125.mat


size(vta_t10_ht05e)
size(vta_t10_ht025e)
size(vta_t10_ht0125e)
size(vta_t10_ht06125e)
E1 = norm(vta_t10_ht05e - vta_t10_ht025e(1:2:end),2)
E2 = norm(vta_t10_ht025e(1:2:end) - vta_t10_ht0125e(1:4:end),2)

conv = (log(E1)-log(E2))/log(2)

EE1 = norm(vta_t10_ht025e - vta_t10_ht0125e(1:2:end),2)
EE2 = norm(vta_t10_ht0125e(1:2:end) - vta_t10_ht06125e(1:4:end),2)

conv2 = (log(EE1)-log(EE2))/log(2)
