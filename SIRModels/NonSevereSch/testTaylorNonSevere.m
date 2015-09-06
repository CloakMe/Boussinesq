clear;

S_0 = 760;
I_0 = 4;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
[ a, r, ro ] = GetParams();

order = 4;
ic = GenerateDerivativesNonSevere( order, IC );

Tend = 180;
t = 0.02;
%Tend = 0.2;
tic
[ T, Y, Sol ] = SIRModelTaylorNonSevere( t, Tend, ic, order );
toc
ResultPlot( 4, T, Y(:,1), '-b', 'Susceptibles' )
ResultPlot( 5, T, Y(:,2), '-r', 'Infected' )
%ResultPlot( 6, T, Y(:,3), '-g', 'Recovered' )

derDeath = a*Y(:,2)';
ResultPlot( 9, T, a*Y(:,2), '-k', 'Recovered Numerical Derivative' )

derDeathAnalytical = DerRemoveAnalitical( T, S_0, S_0 + I_0 + R_0 );
ResultPlot( 10, T, derDeathAnalytical, '-c', 'Recovered Analytical Derivative' )

difff = abs( derDeath - derDeathAnalytical );
ResultPlot( 11, T, difff, '-c', 'Difference' )


