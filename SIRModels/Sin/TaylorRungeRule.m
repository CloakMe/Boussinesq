clear;

u_0 = 0;
v_0 = -1;
IC = [ u_0, v_0 ];

order = 4;
ic = GenerateDerivatives( order, IC );
Tend = 180;
t1 = 0.04;

[ T1, Y1, Sol1 ] = SinModelTaylor( t1, Tend, ic, order );

t2 = 0.02;
[ T2, Y2, Sol2 ] = SinModelTaylor( t2, Tend, ic, order );

t3 = 0.01;
[ T3, Y3, Sol3 ] = SinModelTaylor( t3, Tend, ic, order );

%T3 = 0:t3:Tend;
%da = DerRemoveAnalitical( T3, S_0, S_0 + I_0 + R_0 );

fine1 = Y1(:,1)';
fine2 = Y2(1:2:end,1)';
finer1 = Y2(:,1)';
finer2 = Y3(1:2:end,1)';
%finer2 = da(1:2:end);

fine = (fine1 - fine2);
finer = (finer1 - finer2);

norm0402_L2 = sqrt(t2)*norm( fine(:), 2 );
norm0201_L2 = sqrt(t3)*norm( finer(:), 2 );


fprintf( '||y_04 - y_02||_L2 = %.8e \n', norm0402_L2 );
fprintf( '||y_02 - y_01||_L2 = %.8e \n', norm0201_L2 );
conv_L2 = log( abs( norm0402_L2 / norm0201_L2 ) ) / log(2);
fprintf( 'conv_L2 = %.8e \n', conv_L2 );





