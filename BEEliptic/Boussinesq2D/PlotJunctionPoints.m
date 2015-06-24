function PlotJunctionPoints(x,y,newBigIC)
    quater = floor((length(x)+1)/4);
    figure(2)
    mesh(x(end - quater - 10:1:end - quater + 10),y,newBigIC(end - quater - 10:1:end - quater + 10,:)');
    %mesh(x,y,newBigIC');
    xlabel('x');ylabel('y');
    title('Junction Points');
end