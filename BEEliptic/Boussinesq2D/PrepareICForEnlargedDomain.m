function [bigU,bigUTimeDerivative,P,U,newBigIC,solutionNorms,...
        theta,c1,c2,zeroX,zeroY,tauVector,angl]=...
        PrepareICForEnlargedDomain(oldBigU,compBox,prmtrs,al,bt1,bt2,c,c1,lastTheta,derivative)
    
    x=compBox.x_st2:prmtrs.h:compBox.x_end2; y=compBox.y_st2:prmtrs.h:compBox.y_end2;
    [zeroX,zeroY]=GetZeroNodes(x,y);
    [X,Y]=Domain(x,y);

    bt = bt1/bt2;
    c12 = 1-c^2/bt;
    newBigIC = c1*lastTheta*(c12*X.^2-Y.^2)./(c12*X.^2+Y.^2).^2;
    %newBigIC = c1*lastTheta*(c12^2* X.^4 - 6*c12 * X.^2 .* Y.^2 + Y.^4)./(c12*X.^2+Y.^2).^4;
    points = FindOldGrid(compBox.x_st,compBox.x_end,compBox.y_st,compBox.y_end,x,y);
    newBigIC(points(1):points(2),points(3):points(4)) = oldBigU;

    NewIC = newBigIC(zeroX:end,zeroY:end);
    thet = abs(NewIC(1,1));
    NewIC = NewIC/thet;
    PlotJunctionPoints(x,y,newBigIC);

    [bigU,bigUTimeDerivative,P,U,theta,c1,c2,solutionNorms,tauVector,angl] =...
       sol_ch_v9(NewIC,x,y,prmtrs,bt1,bt2,al,c,thet,zeroX,zeroY,derivative);
       
end