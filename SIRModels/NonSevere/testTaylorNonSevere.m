clear;

S_0 = 9000;
I_0 = 10;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
[ a, r, ro ] = GetParams()

ic = GenerateDerivativesNonSevere( 4, IC );
Tend = 140;
t = 0.02;
%Tend = 0.2;
tic
[ T, Y, Sol ] = SIRModelTaylorNonSevere( t, Tend, ic );
toc
ResultPlot( 4, T, Y(:,1), '-b', 'Susceptibles' )
ResultPlot( 5, T, Y(:,2), '-r', 'Infectives' )
ResultPlot( 6, T, Y(:,3), '-g', 'Death' )

deathAnalytical = RemoveAnalitical( T, S_0, S_0 + I_0 + R_0 );
ResultPlot( 10, T, deathAnalytical, '-c', 'Death Analytical' )

difff = abs( Y(:,3) - deathAnalytical' );
ResultPlot( 11, T, difff, '-c', 'Difference' )