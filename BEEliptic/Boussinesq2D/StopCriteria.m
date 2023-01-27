function [ stopFlag, Px, Py ] = StopCriteria(x, y, zeroX, zeroY, U, axOld, ayOld, minResidual, eps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %eps = 1.000e-010;
    
    midX = floor(length(x(zeroX:end))/2);
    %midY = floor(length(y(zeroY:end))/2);
    
    newX = x(midX+zeroX-1:end);
    %newY = y(midY+zeroY-1:end);
    
    %Py = polyfit(newY, (newY.^2).*U(1,midY:end), 1);
Py=0;
    Px = polyfit(newX', (newX'.^2).*U(midX:end,1), 1);
   
    if( abs(Px(1) - axOld) < eps ) %&& abs(Py(1) - ayOld) < eps 
        stopFlag = 1;
    else
        stopFlag = 0;
    end
    
    if(minResidual > eps)
       stopFlag = 0;
       %Px = [axOld 0];
       %Py = [ayOld 0];
       %return;
    end
    
end

