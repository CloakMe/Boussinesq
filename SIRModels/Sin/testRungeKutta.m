clear;

u_0 = 0;
v_0 = -1;
IC = [ u_0, v_0 ];
Tend = 14000;
t = 0.2;
time_span = 0:t:Tend;
tic
options = odeset( 'RelTol', 1e-4, 'AbsTol', [1e-4 1e-4 1e-5] );
Y = ode4(@SinModel,time_span,IC);
toc
ResultPlot( 1, time_span(end-50:end), Y(end-50:end,1), '-b', 'sin' )
ResultPlot( 2, time_span(end-50:end), Y(end-50:end,2), '-r', '-cos' )
