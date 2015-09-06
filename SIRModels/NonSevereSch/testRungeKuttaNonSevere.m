clear;

S_0 = 760;
I_0 = 4;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
[ a, r, ro ] = GetParams();

Tend = 180;
t = 0.02;
time_span = 0:t:Tend;
T = time_span;
tic
options = odeset( 'RelTol', 1e-4, 'AbsTol', [1e-4 1e-4 1e-5] );
Y = ode4(@SIRModelNonSevere,T,IC);
toc
ResultPlot( 1, T, Y(:,1), '-b', 'Susceptibles' )
ResultPlot( 2, T, Y(:,2), '-r', 'Infectives' )
ResultPlot( 3, time_span, Y(:,3), '-g', 'Recovered' )


ResultPlot( 9, T, a * Y(:,2), '-k', 'Recovered Numerical Derivative' )

derDeathAnalytical = DerRemoveAnalitical( T, S_0, S_0 + I_0 + R_0 );
ResultPlot( 10, T, derDeathAnalytical, '-c', 'Recovered Analytical Derivative' )

difff = abs( a * Y(:,2) - derDeathAnalytical' );
ResultPlot( 11, T, difff, '-c', 'Difference' )


