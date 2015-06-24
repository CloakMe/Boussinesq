clear;
h = .1;
st = 15;
x = -st:h:st;
y = x;   
sx = size(x, 2)


[ ff, dhb ] = BEUtilities.GetFinDiffMat( sx, h );
domUtils = BEDomainUtils( x, y, 4 );
M = rand( sx,sx );
fd = BEUtilities.GetFinDiffCoeff( [-2,-1,0,1,2], 2 );
tic
[ result1 ] = domUtils.DeltaH( M, fd' );
toc
tic
result2 = dhb*M + M*dhb;
toc