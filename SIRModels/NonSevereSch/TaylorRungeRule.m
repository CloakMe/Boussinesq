clear;

S_0 = 760;
I_0 = 4;
R_0 = 0;
IC = [ S_0, I_0, R_0 ];
[ a, r, ro ] = GetParams();

order = 4;
ic = GenerateDerivativesNonSevere( order, IC );

Tend = 180;
t1 = 0.04;
[ T1, Y1, Sol1 ] = SIRModelTaylorNonSevere( t1, Tend, ic, order );

t2 = 0.02;
[ T2, Y2, Sol2 ] = SIRModelTaylorNonSevere( t2, Tend, ic, order );

t3 = 0.01;
[ T3, Y3, Sol3 ] = SIRModelTaylorNonSevere( t3, Tend, ic, order );

%T3 = 0:t3:Tend;
%da = DerRemoveAnalitical( T3, S_0, S_0 + I_0 + R_0 );

fine1 = a * Y1(:,2)';
fine2 = a * Y2(1:2:end,2)';
finer1 = a * Y2(:,2)';
finer2 = a * Y3(1:2:end,2)';
%finer2 = da(1:2:end);

fine = (fine1 - fine2);
finer = (finer1 - finer2);

norm0402_L2 = sqrt(t2)*norm( fine(:), 2 );
norm0201_L2 = sqrt(t3)*norm( finer(:), 2 );


fprintf( '||y_04 - y_02||_L2 = %.8e \n', norm0402_L2 );
fprintf( '||y_02 - y_01||_L2 = %.8e \n', norm0201_L2 );
conv_L2 = log( abs( norm0402_L2 / norm0201_L2 ) ) / log(2);
fprintf( 'conv_L2 = %.8e \n', conv_L2 );

derDeath1 = a*Y1(:,2);
ResultPlot( 12, T1, derDeath1, '-k', 'Recovered Numerical Derivative' )

derDeath2 = a*Y2(:,2);
ResultPlot( 13, T2, derDeath2, '-k', 'Recovered Numerical Derivative' )

derDeath3 = a*Y3(:,2);
ResultPlot( 14, T3, derDeath3, '-k', 'Recovered Numerical Derivative' )



