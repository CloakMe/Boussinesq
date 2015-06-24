clear;
h = .1;
x = -3:h:3;
y = -3:h:3;
order = 2; 

domUtils = BEDomainUtils( x, y, order );

[X,Y] = domUtils.GetNet(x,y);
v = BEHeatMassTransferEquation( X,Y );
mesh( x, y, v' );