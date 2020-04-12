function [bigU,bigUTimeDerivative,P,U,thetaCont,c1,c2,solutionNormsCont,tauVecCont,anglCont,nx,ny,h] =...
    PreSolverFromExistingSol(x,y,U,compBox,prmtrs,bt1,bt2,al,c,th,derivative,h)

    tic
    prmtrs.h = h;
    prmtrs.tau = prmtrs.tau/10;
    nx = compBox.x_st2:h:compBox.x_end2;
    ny = compBox.y_st2:h:compBox.y_end2; 
    [zeroX,zeroY]=GetZeroNodes(x,y);
    [X,Y]=Domain(x(zeroX:end),y(zeroY:end));
    [zeronX,zeronY]=GetZeroNodes(nx,ny);
    [nX,nY]=Domain(nx(zeronX:end),ny(zeronY:end));
  
    nU = griddata(X,Y,U,nX,nY,'cubic');
 
    
    fprintf('elapsed pre solver time = %d \n', toc);
    fprintf('new solution (x,y) size = (%d,%d)\n', length(nx), length(ny));
    
    tic
    [bigU,bigUTimeDerivative,P,U,thetaCont,c1,c2,solutionNormsCont,tauVecCont,anglCont] =...
       sol_ch_v8(nU,nx,ny,prmtrs,bt1,bt2,al,c,th,zeronX,zeronY,derivative);   
    fprintf('elapsed solver time = %d \n', toc);
end