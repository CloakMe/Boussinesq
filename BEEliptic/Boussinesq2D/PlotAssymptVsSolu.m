function PlotAssymptVsSolu( x, y, h, bigU, muTheta, c)

    shift = ceil(2/h);
    [zeroX,zeroY]=GetZeroNodes(x,y);
    stX = zeroX+shift-1;
    stY = zeroY+shift-1;
    newX = x(stX:end);
    newY = y(stY:end);
    
    assymptYeqZero = + muTheta ./ ( (1-c^2) * newX.^2 ); % 
    assymptXeqZero = - muTheta ./ newY.^2; % 

    figure(5)
    plot(newY, (newY.^2).*bigU(zeroX,stY:end), 'b', newY, assymptXeqZero .* newY.^2, 'k' ) %(newY.^2).*
    xlabel('y')
    title('x==0 Cross section * y^2');
    
    figure(6)
    plot(newX, (newX'.^2).*bigU(stX:end,zeroY), 'b', newX, assymptYeqZero .* newX.^2, 'k' ) % (newX'.^2).*
    xlabel('x')
    title('y==0 Cross section * x^2');
    
    figure(7)
    plot(newY, bigU(zeroX,stY:end), 'b', newY, assymptXeqZero, 'k' ) %(newY.^2).*
    xlabel('y')
    title('x==0 Cross section');
    
    figure(8)
    plot(newX, bigU(stX:end,zeroY), 'b', newX, assymptYeqZero, 'k' ) % (newX'.^2).*
    xlabel('x')
    title('y==0 Cross section');
end