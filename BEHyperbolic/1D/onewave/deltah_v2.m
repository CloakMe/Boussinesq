function [resmat]=deltah_v2(sx,h)

resmat=-2*diag(ones(sx,1));
resmat(1,1)=1;
resmat(1,2)=-2;
resmat(1,3)=1;
for l=2:sx-1
    resmat(l,l-1)=1;
    resmat(l,l+1)=1;
end
resmat(sx,sx-2)=1;
resmat(sx,sx-1) =-2;
resmat(sx,sx) =1;
%resmat=resmat/h^2;