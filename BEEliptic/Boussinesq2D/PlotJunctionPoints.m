function PlotJunctionPoints(x,y,newBigIC,quater)
    
    figure(2)
    mesh(x(end - quater - 10:1:end - quater + 10),y,newBigIC(end - quater - 10:1:end - quater + 10,:)');
    %mesh(x,y,newBigIC');
    xlabel('x');ylabel('y');
    title('Junction Points');
end