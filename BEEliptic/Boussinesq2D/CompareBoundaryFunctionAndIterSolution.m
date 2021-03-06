function [uBvsU_L2,PBvsP_L2]=CompareBoundaryFunctionAndIterSolution(U,P,approximationBoundaryU,approximationBoundaryP,step,h)
    [xLength,yLength] = size(U);
    [yLen,xLen,yHei,xHei] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

    ue_i=[U(xHei,yLen) U(xLen,yHei)'];  
    pe_i=[P(xHei,yLen) P(xLen,yHei)'];
    
    uBvsU_L2 =  h*norm( approximationBoundaryU(:) - ue_i(:), 2);
    PBvsP_L2 =  h*norm( approximationBoundaryP(:) - pe_i(:), 2);
end