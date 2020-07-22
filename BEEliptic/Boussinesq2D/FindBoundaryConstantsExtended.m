function [muU1,muU2,muP1,muP2]=FindBoundaryConstantsExtended(U,P, ...
    approximationBoundaryU1, ...
    approximationBoundaryU2, ...
    approximationBoundaryP1, ...
    approximationBoundaryP2, ...
    step)
    [xLength,yLength] = size(U);
    [yLen,xLen,yHei,xHei] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

    ue_i=[U(xHei,yLen) U(xLen,yHei)'];  
    pe_i=[P(xHei,yLen) P(xLen,yHei)'];
    
    al_1_1 = sum(sum(approximationBoundaryU1.*approximationBoundaryU1));
    al_1_2 = sum(sum(approximationBoundaryU1.*approximationBoundaryU2));
    b_1 = sum(sum(ue_i.*approximationBoundaryU1));
    
    al_2_1 = al_1_2;
    al_2_2 = sum(sum(approximationBoundaryU2.*approximationBoundaryU2));    
    b_2 = sum(sum(ue_i.*approximationBoundaryU2));
    A = [al_1_1 al_1_2; al_2_1 al_2_2];
    B = [b_1 b_2];
    
    result = linsolve(A, B');
    muU1 = result(1);
    muU2 = result(2);

    if(muU2 < 0)
        muU2 = 0;
    end
    
    al_1_1 = sum(sum(approximationBoundaryP1.*approximationBoundaryP1));
    al_1_2 = sum(sum(approximationBoundaryP1.*approximationBoundaryP2));
    b_1 = sum(sum(pe_i.*approximationBoundaryP1));
    
    al_2_1 = al_1_2;
    al_2_2 = sum(sum(approximationBoundaryP2.*approximationBoundaryP2));    
    b_2 = sum(sum(pe_i.*approximationBoundaryP2));
    A = [al_1_1 al_1_2; al_2_1 al_2_2];
    B = [b_1 b_2];
    
    result = linsolve(A, B');
    muP1 = result(1);
    muP2 = result(2);

    if(muP2 < 0)
        muP2 = 0;
    end
end