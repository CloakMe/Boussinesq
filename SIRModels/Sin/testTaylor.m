clear;

u_0 = 0;
v_0 = -1;
IC = [ u_0, v_0 ];
Tend = 14000;
t = 0.2;

order = 4;
ic = GenerateDerivatives( order, IC );

%Tend = 0.2;
tic
[ T, Y, Sol ] = SinModelTaylor( t, Tend, ic, order );
toc

ResultPlot( 4, T(end-50:end), Y(end-50:end,1), '-b', 'sin taylor' )
ResultPlot( 5, T(end-50:end), Y(end-50:end,2), '-r', '-cos taylor' )