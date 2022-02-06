clear;clc;     
ff = load('dtvta_t10_h025_t025.dat'); 
fv = load('vta_t10_h025_t025.dat'); 
     %size(ff)
     
     load vta_t10_h025_t025.mat
     -10.01 + tt(end)
     size(tt)

        vz = fv(:,end);
        sx = size(x,2)*size(x,1);
        dvz = ff(:,end);
        plot(x,dvz,x,vz)
        h = x(2) - x(1);
[d2vz, d3vz, d4vz, d5vz] = calc_der(vz,dvz,[1 -2 1],[-1 (h^2 + 2) -1],h,sx,3,3);
        time= 0.01151011402222;
        vta_t10_h025_t025e  =  vz +  time*dvz + time^2*d2vz/2 + time^3*d3vz/6 + time^4*d4vz/24;

        save vta_t10_h025_t025.mat vta_t10_h025_t025e tt x;