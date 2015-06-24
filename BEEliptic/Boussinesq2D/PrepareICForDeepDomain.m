function [bigU,bigUTimeDerivative,P,U,newBigIC,solutionNorms,...
        theta,c1,c2,zeroX,zeroY,tauVector,angl]=...
        PrepareICForDeepDomain(oldBigU,compBox,prmtrs,al,bt1,bt2,c,c1,lastTheta,derivative)
    %DEEP
    finer = 2;
    xNew=compBox.x_st:prmtrs.h/finer:compBox.x_end;
    yNew=compBox.y_st:prmtrs.h/finer:compBox.y_end;
    
x2 = x(1):h/finer:x(end);
y2 = y(1):h/finer:y(end);
x2(1:finer:end) = x;
y2(1:finer:end) = y;
hn = h/finer;
 ox1 = ones(1,sy);
    X = x'*ox1;
   if(sx~=sy)
        oy1 = ones(1,sx);
        Y = (y'*oy1)';
    else
        Y = (y'*ox1)';
   end
   
 ox1 = ones(1,length(y2));
    X2 = x2'*ox1;
   if(sx~=sy)
        oy1 = ones(1,length(x2));
        Y2 = (y2'*oy1)';
   else
        Y2 = (y2'*ox1)';
   end
   
[y3,x3,ZI] = griddata(X,Y,bU,X2,Y2,'linear');
    
    [zeroX,zeroY]=GetZeroNodes(x,y);
    [X,Y]=Domain(x,y);

    c12 = 1-c^2;
    newBigIC = c1*lastTheta*(c12*X.^2-Y.^2)./(c12*X.^2+Y.^2).^2;
    %newBigIC = c1*lastTheta*(c12^2* X.^4 - 6*c12 * X.^2 .* Y.^2 + Y.^4)./(c12*X.^2+Y.^2).^4;
    points = FindOldGrid(compBox.x_st,compBox.x_end,compBox.y_st,compBox.y_end,x,y);
    newBigIC(points(1):points(2),points(3):points(4)) = oldBigU;

    NewIC = newBigIC(zeroX:end,zeroY:end);
    thet = abs(NewIC(1,1));
    NewIC = NewIC/thet;
    PlotJunctionPoints(x,y,newBigIC);

    [bigU,bigUTimeDerivative,P,U,theta,c1,c2,solutionNorms,tauVector,angl] =...
       sol_ch_v8(NewIC,x,y,prmtrs,bt1,bt2,al,c,thet,zeroX,zeroY,derivative);
       
end