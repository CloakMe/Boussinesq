function [sw ]= IsAlgoDiverging(subCounter,residualInfNorm)
    sw = 0;

    if( residualInfNorm(subCounter)>40000 || ...
        isnan(residualInfNorm(subCounter)) || ...
        isinf(residualInfNorm(subCounter)))

            warning(' Algo diverges!');
            sw = 1;
    end
end