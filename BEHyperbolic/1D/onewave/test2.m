clear;clc;
load vpdt_t10_h1.mat
load vpdt_t10_h05.mat
load vpdt_t10_h025.mat
load vpdt_t10_h0125.mat

%load sv0025.mat
%load sv00125.mat
%load sv000625.mat


size(vpdt_t10_h1e)
size(vpdt_t10_h05e)
size(vpdt_t10_h025e)
size(vpdt_t10_h0125e)
E1 = norm(vpdt_t10_h1e - vpdt_t10_h05e(1:2:end),2)
E2 = norm(vpdt_t10_h05e(1:2:end) - vpdt_t10_h025e(1:4:end),2)

conv = (log(E1)-log(E2))/log(2)

EE1 = norm(vpdt_t10_h05e - vpdt_t10_h025e(1:2:end),2)
EE2 = norm(vpdt_t10_h025e(1:2:end) - vpdt_t10_h0125e(1:4:end),2)

conv2 = (log(EE1)-log(EE2))/log(2)

