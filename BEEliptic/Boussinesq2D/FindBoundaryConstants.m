function [c1,c2]=FindBoundaryConstants(U,P,approximationBoundaryU,approximationBoundaryP,step)
    [xLength,yLength] = size(U);
    [yLen,xLen,yHei,xHei] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

    ue_i=[U(xHei,yLen) U(xLen,yHei)'];  
    pe_i=[P(xHei,yLen) P(xLen,yHei)'];
    
    c1 =  sum(sum(ue_i.*approximationBoundaryU)) /sum(sum(approximationBoundaryU.^2));
    if(c1 < 0)
        c1 = 0;
    end
    c2 =  sum(sum(pe_i.*approximationBoundaryP)) /sum(sum(approximationBoundaryP.^2));
    if(c2 < 0)
        c2 = 0;
    end
end