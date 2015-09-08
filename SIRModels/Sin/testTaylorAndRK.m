clear;

u_0 = 0;
v_0 = -1;
IC = [ u_0, v_0 ];
Tend = 14000;
t = 0.1;

time_span = 0:t:Tend;
tic
Yr = ode4(@SinModel,time_span,IC);
toc

order = 4;
ic = GenerateDerivatives( order, IC );

%Tend = 0.2;
tic
[ T, Yt, Sol ] = SinModelTaylor( t, Tend, ic, order );
toc




diff_r = Yr(end-50:end,1) + sin( T(end-50:end) )';
ResultPlot( 1, T(end-50:end), abs( diff_r ) , '-k', 'diff runge' );

diff_t = Yt(end-50:end,1) + sin( T(end-50:end) )';
ResultPlot( 4, T(end-50:end), abs(diff_t), '-b', 'diff taylor' )
return;
diff_r1 = Yr(:,1) + sin( T )';
ResultPlot( 2, T, abs( diff_r1 ) , '-k', 'diff runge' );

diff_t1 = Yt(:,1) + sin( T )';
ResultPlot( 5, T, abs(diff_t1), '-b', 'diff taylor' )