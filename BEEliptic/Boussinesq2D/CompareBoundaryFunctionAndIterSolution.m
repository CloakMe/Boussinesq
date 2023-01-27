function [uBvsU_L2, PBvsP_L2, U_L2, P_L2, uBvsU_Inf, PBvsP_Inf, U_Inf, P_Inf] = ...
    CompareBoundaryFunctionAndIterSolution(U, P, approximationBoundaryU, approximationBoundaryP, step, h)
    [xLength,yLength] = size(U);
    [yLen,xLen,yHei,xHei] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

    ue_i=[U(xHei,yLen) U(xLen,yHei)'];  
    pe_i=[P(xHei,yLen) P(xLen,yHei)'];
    try
        uBvsU_L2 =  h*norm( approximationBoundaryU(:) - ue_i(:), 2);
        PBvsP_L2 =  h*norm( approximationBoundaryP(:) - pe_i(:), 2);
        U_L2 =  h*norm( ue_i(:), 2);
        P_L2 =  h*norm( pe_i(:), 2);
    catch
        U_L2 =  0;
        P_L2 =  0;
    end
    
    uBvsU_Inf = max( abs( approximationBoundaryU(:) - ue_i(:) ) );
    PBvsP_Inf = max( abs( approximationBoundaryP(:) - pe_i(:) ) );
    U_Inf = max( abs( ue_i(:) ) );
    P_Inf = max( abs( pe_i(:) ) );
end