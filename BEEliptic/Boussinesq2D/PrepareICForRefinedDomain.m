function [bigU,bigUTimeDerivative,P,U,newBigU,solutionNorms,...
        theta,c1,c2,zeroX,zeroY,tauVector,angl]=...
        PrepareICForRefinedDomain(oldBigU,compBox,prmtrs,al,bt1,bt2,c,c1,lastTheta,derivative)

    finer = 2;
    x=compBox.x_st:prmtrs.h:compBox.x_end;
    y=compBox.y_st:prmtrs.h:compBox.y_end;
    xNew=compBox.x_st:prmtrs.h/finer:compBox.x_end;
    yNew=compBox.y_st:prmtrs.h/finer:compBox.y_end;
    
    [X,Y]=Domain(x,y);
    [XNew,YNew]=Domain(xNew,yNew);
    
    [y3,x3,newBigU] = griddata(X,Y,oldBigU,XNew,YNew,'cubic');
    PlotRefinedDifference(xNew,yNew,finer,newBigU,oldBigU);
    
    [zeroX,zeroY]=GetZeroNodes(xNew,yNew);
    NewIC = newBigU(zeroX:end,zeroY:end);
    thet = abs(NewIC(1,1));
    NewIC = NewIC/thet;
    
    prmtrs.h = prmtrs.h/finer;
    [bigU,bigUTimeDerivative,P,U,theta,c1,c2,solutionNorms,tauVector,angl] =...
       sol_ch_v8(NewIC,xNew,yNew,prmtrs,bt1,bt2,al,c,thet,zeroX,zeroY,derivative);
end