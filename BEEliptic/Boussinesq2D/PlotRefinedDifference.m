function PlotRefinedDifference(x,y,finer,newBigIC,oldBigU)    

    figure(5)
    mesh(x,y,newBigIC(1:end,1:end)');%-oldBigU'
    %mesh(x,y,newBigIC');
    xlabel('x');ylabel('y');
    title('Difference b/n refined and non-refined solution');
end