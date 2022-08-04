clear;
addpath('..\..\..\Common');
h = .1;
st = 30.0;
x = -st+h:h:st-h;
y = x;   
sx = size(x, 2)

[X,Y]=Domain(x,y);

[ dhb ] = BEUtilities.GetFinDiffMatZeroBnd( sx, 6 );
domUtils = BEDomainUtils( x, y, 6 );
M = exp(- (X.^2 + Y.^2) );
figure(1)
mesh(x,x,(M)')

fd = BEUtilities.GetFinDiffCoeff( [-3,-2,-1,0,1,2,3], 2 )';
zeroMatrix = zeros(size(M));

fprintf('DeltaH: \n');
tic
[ result01 ] = domUtils.DeltaH( M, fd )/h^2;
toc
%fastest
fprintf('DeltaHZeroBnd: \n');
tic
[ result02 ] = domUtils.DeltaHZeroBnd( M )/h^2;
toc
%max(max(M))
fprintf('diff %f: \n', max(max(result01-result02)));

%fastest
fprintf('DeltaHEvenFunZeroBnd: \n');
tic
[ result11 ] = domUtils.DeltaHEvenFunZeroBnd(M)/h^2;
toc
fprintf('diff %f: \n', max(max(result01-result11)));


fprintf('DeltaEvenFunctions: \n');
tic
[ result12 ] = YDerivativeEvenFunctions(M,zeroMatrix,0,fd) + XDerivativeEvenFunctions(M,zeroMatrix,0,fd);
toc
fprintf('diff %f: \n', max(max(result12-result11)));

fprintf('\nfull matrix: \n');
tic
result3 = (dhb*M + M*dhb')/h^2;
toc
fprintf('diff %f: \n', max(max(result3-result02)));
return;
result4 = (A*M + M*A')/h^2;
fprintf('diff %f: \n', max(max(result3-result4)));


figure(2)
mesh(x,x,(result3-result4)')

figure(3)
mesh(x,x,(result02-result4)')




fprintf('\nBandMatMult(*,*,*): \n');
tic
res3 = BEUtilities.BandMatMult(fd,M,12) + BEUtilities.BandMatMult(fd,M,12);
toc

fprintf('\nBandMatMult(*,*): \n');
tic
res4 = BEUtilities.BandMatMult(fd,M) + BEUtilities.BandMatMult(fd,M);
toc

[ dmat ] = BEUtilities.GetFinDiffMatZeroBnd( 10, 6 );
[S,D] = eig(-dmat)
[invS, invD]=eig(inv(-dmat))