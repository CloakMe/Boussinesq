clear;
n=10;

%h = 1/n;
h = this.h;
n = this.sx;
W = zeros(n,n);
L = zeros(n,1);
for i=1:n
    L(i) = 4*sin( i*pi/80*h/2)^2;
    for j=1:n
        W(i,j) = sin(i*j*pi/80*h);
    end
end

figure(1)
plot( 1:this.sx, L, 'r', 1:this.sx, diag(w), 'k')

figure(2)
plot(1:this.sx, L-diag(w))

figure(3)
mesh(1:this.sx, 1:this.sx, this.eigenFinDiffMat)

figure(4)
mesh(1:this.sx, 1:this.sx, W - this.eigenFinDiffMat)

figure(5)
plot(1:this.sx, W(1:end,50)/norm(W(1:end,50),2), 1:this.sx,this.eigenFinDiffMat(1:end,50), 'r')

figure(6)
plot(1:this.sx, W(1:end,50)/norm(W(1:end,50),2) - this.eigenFinDiffMat(1:end,50))


dhb = BEUtilities.GetFinDiffMat( 2*n+1, 2, h );
F = (-dhb)*W' - W'*diag(L)
%A=W'*diag(L)*W
    