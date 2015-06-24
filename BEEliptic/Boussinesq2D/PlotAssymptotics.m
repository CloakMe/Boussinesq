function PlotAssymptotics(x,y,h,zeroX,zeroY,bigU,holdon)
    if(nargin == 6)
        holdon = 0;
    end
    shift = ceil(6/h);
    stX = zeroX+shift-1;
    stY = zeroY+shift-1;
    newX = x(stX:end);
    newY = y(stY:end);
if(size(x,2) == size(bigU,1))
    %{
    figure(1)
    loglog(y(zeroY:end), abs( bigU(zeroX,zeroY:end) ) ) %(1+end)/2
    xlabel('y')
    title('x==0 Cross section Log-Log');
    figure(2)
    loglog(x(zeroX:end),abs( bigU(zeroX:end,zeroY) ) )
    xlabel('x')
    title('y==0 Cross section Log-Log');
    %}
  
    figure(3)
    HoldOn(holdon);
    plot(newY, (newY.^2).* bigU(zeroY,stY:end) ) %(1+end)/2
    HoldOff(holdon);
    xlabel('y')
    title('x==0 Cross section');
    
    figure(4)
    HoldOn(holdon);
    plot(newX,(newX'.^2).*bigU(stX:end,zeroY) )
    HoldOff(holdon);
    xlabel('x')
    title('y==0 Cross section');
else
    figure(3)
    HoldOn(holdon);
    plot(newY, (newY.^2).* bigU(1,shift:end) ) %(1+end)/2
    HoldOff(holdon);
    xlabel('y')
    title('x==0 Cross section');
    figure(4)
    HoldOn(holdon);
    plot(newX,(newX'.^2).*bigU(shift:end,1) )
    HoldOff(holdon);
    xlabel('x')
    title('y==0 Cross section');
end
    
end

function HoldOn(holdon)
    if(holdon ~= 0)
        hold on;
    end
end

function HoldOff(holdon)
    if(holdon ~= 0)
        hold off;
    end
end