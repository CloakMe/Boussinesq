function PlotBoundary(x,y,zeroX, U, outerTopBoundaryF)
    
    augBoundarySize = 0;
    if(nargin == 4)
        BND = U(:,end-7:end);
        augBoundarySize = 4;
    else
        BND = [U(:,end-7:end) outerTopBoundaryF];
    end
    
    hy = y(2)-y(1);
    figure(4);
    mesh(x(zeroX:end),y(end-11+augBoundarySize:end)+(4-augBoundarySize)*hy,BND');
    xlabel('x');    ylabel('y');
    view(-13,16);
    title('Boundary Smoothness');
end
    %x2 = [x x(end)+h x(end)+2*h x(end)+3*h x(end)+4*h];
    %y2 = [y y(end)+h y(end)+2*h y(end)+3*h y(end)+4*h];
    
    %mesh(x(zeroX:end),y(zeroX:end),U');