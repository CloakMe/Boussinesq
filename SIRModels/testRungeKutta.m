clear;

S_0 = 762;
I_0 = 1;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
tic
options = odeset( 'RelTol', 1e-4, 'AbsTol', [1e-4 1e-4 1e-5] );
[T,Y] = ode45(@SIRModel,[0, 14],IC,options);
toc
ResultPlot( 1, T, Y(:,1), '-y', 'Susceptibles' )
ResultPlot( 2, T, Y(:,2), '-r', 'Infectives' )
ResultPlot( 3, T, Y(:,3), '-g', 'Recovered' )
