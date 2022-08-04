function [W]=eigvect(sx)
h=1/(sx+1);
i=1:sx;
j=i;
W = sqrt(2*h)*sin((pi*h)*i'*j);