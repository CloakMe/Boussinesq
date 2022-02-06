function [resmat]=deltah(sx)

resmat=-2*diag(ones(sx,1));
resmat(1,2)=1;
for l=2:sx-1
    resmat(l,l-1)=1;
    resmat(l,l+1)=1;
end
resmat(sx,sx-1)=1;
%resmat=resmat/h^2;