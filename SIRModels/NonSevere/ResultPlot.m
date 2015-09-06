function ResultPlot( picNum, x , y, lineType, ResultTitle )
    figure( picNum );
    plot( x, y, lineType );
    xlabel( 't' ); ylabel( 'Y' ); 
    title( ResultTitle );
end

