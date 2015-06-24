function PlotResidual(x,y,crrntResidual)
    figure(3)
    mesh(x,y,crrntResidual');
    xlabel('x');    ylabel('y');
    view(-13,30);
    title('Residual');
end
%{
    figure(7)
    mesh(x(zeroX:end),y(zeroY:end),U');
    xlabel('x');    ylabel('y');
    view(-13,30);
    title('current Sol');
%}