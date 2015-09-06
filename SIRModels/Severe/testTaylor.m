clear;

S_0 = 762;
I_0 = 1;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
ic = GenerateDerivatives( 4, IC );
Tend = 14;
t = 0.05;
%Tend = 0.2;
tic
[ T, Y, Sol ] = SIRModelTaylor( t, Tend, ic );
toc
ResultPlot( 4, T, Y(:,1), '-b', 'Susceptibles' )
ResultPlot( 5, T, Y(:,2), '-r', 'Infectives' )
ResultPlot( 6, T, Y(:,3), '-g', 'Recovered' )