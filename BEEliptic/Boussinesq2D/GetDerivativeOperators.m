function [derivative,order]=GetDerivativeOperators(pointsLocation,h)

    firstDerivative = GetFiniteDifferenceCoeff(pointsLocation,1)'/h;
    secondDerivative = GetFiniteDifferenceCoeff(pointsLocation,2)'/h^2;
    derivative = struct('first',{firstDerivative},'second',{secondDerivative});
    order = (length(secondDerivative)-1);
end