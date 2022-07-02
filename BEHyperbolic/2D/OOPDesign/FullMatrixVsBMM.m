clear;
addpath('..\..\..\Common');
h = .1;
st = 10;
x = -st:h:st;
y = x;   
sx = size(x, 2)


[ dhb ] = BEUtilities.GetFinDiffMat( sx, 6 );
domUtils = BEDomainUtils( x, y, 6 );
M = rand( sx,sx );
fd = BEUtilities.GetFinDiffCoeff( [-3,-2,-1,0,1,2,3], 2 )';
zeroMatrix = zeros(size(M));

fprintf('DeltaH: \n');
tic
[ result01 ] = domUtils.DeltaH( M, fd );
toc
%fastest
fprintf('DeltaHZeroBnd: \n');
tic
[ result02 ] = domUtils.DeltaHZeroBnd( M );
toc
max(max(M))
max(max(result01(3:end-2,3:end-2)-result02(3:end-2,3:end-2)))
%fastest
fprintf('DeltaHEvenFunZeroBnd: \n');
tic
[ result11 ] = domUtils.DeltaHEvenFunZeroBnd(M);
toc


fprintf('DeltaEvenFunctions: \n');
tic
[ result12 ] = YDerivativeEvenFunctions(M,zeroMatrix,0,fd) + XDerivativeEvenFunctions(M,zeroMatrix,0,fd);
toc
max(max(M))
max(max(result11-result12))
fprintf('\nfull matrix: \n');
tic
result2 = dhb*M + M*dhb;
toc

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