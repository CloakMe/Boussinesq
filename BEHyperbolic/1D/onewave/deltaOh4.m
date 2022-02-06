function [resmat]=deltaOh4(sx)

resmat=-30*diag(ones(sx,1));
resmat(1,1) = -29;
resmat(sx,sx) = -29;
resmat(1,2)=16;
resmat(1,3)=-1;
resmat(2,4)=-1;
for l=2:sx-1
    resmat(l,l-1)=16;
    resmat(l,l+1)=16;
end
for l=3:sx-2
    resmat(l,l-2)=-1;
    resmat(l,l+2)=-1;
end
resmat(sx,sx-1)=16;
resmat(sx,sx-2)=-1;
resmat(sx-1,sx-3)=-1;
resmat=resmat/12;