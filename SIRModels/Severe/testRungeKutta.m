clear;

S_0 = 762;
I_0 = 1;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
Tend = 14;
t = 0.05;
time_span = 0:t:Tend;
tic
options = odeset( 'RelTol', 1e-4, 'AbsTol', [1e-4 1e-4 1e-5] );
Y = ode4(@SIRModel,time_span,IC);
toc
ResultPlot( 1, time_span, Y(:,1), '-b', 'Susceptibles' )
ResultPlot( 2, time_span, Y(:,2), '-r', 'Infected' )
ResultPlot( 3, time_span, Y(:,3), '-g', 'Recovered' )
