function DrawCrossSections(x,y,h,zeroX,zeroY,al,bt,c,theta,bigU,bigUTimeDer,bigIC,U,compBox,secondDerivative)
  
    fig5=figure(5);
    set(fig5, 'OuterPosition', [0.0      	30.0        380.0     340.0]);
    plot(y(1:zeroY),bigU(zeroX,1:zeroY)); %(1+end)/2
    xlabel('y');
    title('x==0 cross-section');
    %axis([y_st 0 -0.194 -.0]);

    fig6=figure(6);
    set(fig6, 'OuterPosition', [0.0      	30.0        380.0     340.0]);
    plot(x(zeroX:end),bigU(zeroX:end,zeroY));
    xlabel('x');
    title('y==0 cross-section');

    figure(7)
    plot(y,bigU(1,:)) %(1+end)/2
    xlabel('y')
    title('x=x_s_t cross-section');
    figure(8)
    plot(x,bigU(:,1))
    xlabel('x')
    title('y=y_s_t cross-section');

end