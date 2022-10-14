function PlotAssymptVsSolu( x, y, h, bigU, muTheta, c, factor)

    if( nargin == 6 )
        factor = 2;
    end

    shift = ceil(10/h);
    [zeroX,zeroY]=GetZeroNodes(x,y);
    stX = zeroX+shift-1; %zeroX+shift-1;  length(x) - shift;
    stY = zeroY+shift-1; %zeroY+shift-1; length(y) - shift;
    newX = x(stX:end);
    newY = y(stY:end);
        
    assymptYeqZero = + muTheta ./ ( (1-c^2) * newX.^factor ); % 
    assymptXeqZero = - muTheta ./ newY.^factor; % 
    
    
    figure(5)
    plot(newY, (newY.^factor).*bigU(zeroX,stY:end), 'b', newY, assymptXeqZero .* newY.^factor, 'k' ) %(newY.^2).*
    xlabel('y')
    title('x==0 Cross section * y^2');
    
    figure(6)
    plot(newX, (newX'.^factor).*bigU(stX:end,zeroY), 'b', newX, assymptYeqZero .* newX.^factor, 'k' ) % (newX'.^2).*
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