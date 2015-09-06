clear;

S_0 = 9000;
I_0 = 10;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
[ a, r, ro ] = GetParams()


Tend = 140;
t = 0.02;
time_span = 0:t:Tend;
T = time_span;
tic
Y = ode4(@SIRModelNonSevere,T,IC);
toc
ResultPlot( 1, T, Y(:,1), '-b', 'Susceptibles' )
ResultPlot( 2, T, Y(:,2), '-r', 'Infectives' )
ResultPlot( 3, T, Y(:,3), '-g', 'Deaths' )

deathAnalytical = RemoveAnalitical( T, S_0, S_0 + I_0 + R_0 );
ResultPlot( 10, T, deathAnalytical, '-c', 'Death Analytical' )

difff = abs( Y(:,3) - deathAnalytical' );
ResultPlot( 11, T, difff, '-c', 'Difference' )



