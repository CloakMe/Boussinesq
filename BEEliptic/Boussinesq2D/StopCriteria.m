function [ flag, ax, ay ] = StopCriteria(x, y, h, zeroX, zeroY, U, axOld, ayOld, eps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %eps = 1.000e-010;
    midX = floor(length(x(zeroX:end))/2);
    midY = floor(length(y(zeroY:end))/2);
    
    newX = x(midX+zeroX-1:end);
    newY = y(midY+zeroY-1:end);
    
    Py = polyfit(newY, (newY.^2).*U(1,midY:end), 1);

    Px = polyfit(newX', (newX'.^2).*U(midX:end,1), 1);

    ax = Px(1);
    ay = Py(1);
    
    if(abs(ax - axOld) < eps && abs(ay - ayOld) < eps)
        flag = 1;
    else
        flag = 0;
    end
    
end

