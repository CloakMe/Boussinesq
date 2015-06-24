function angl =  Deviation(residualInfNorm,subCounter)
    if(subCounter - 2 <1)
        angl = 0;
        return;
    end
    
    exponent = 10^(-ceil(log10(residualInfNorm(subCounter)))+1);
    p1 = residualInfNorm(subCounter-2);
    p2 = residualInfNorm(subCounter-1);
    p3 = residualInfNorm(subCounter) ;
    
    side1 = sqrt((p2 - p1)^2 + 1/exponent);
    side2 = sqrt((p3 - p2)^2 + 1/exponent);
    
    side1Suare = (p2 - p1)^2 + 1/exponent;
    side2Suare = (p3 - p2)^2 + 1/exponent;
    side3Suare = (p3 - p1)^2 + 4/exponent;
    cos = -(side3Suare - side2Suare - side1Suare)/(2*side1*side2);
    if(p1 + p3 >= 2*p2)
        angl = acos( cos )/pi*180;
    else
        angl = 360 - acos( cos )/pi*180;
    end
    if(~isreal(angl) )
        if(abs(1 + cos ) < 0.0001) % cos =~= -1
            angl = 180;
            return;
        end
        error('angle is complex!');
    end
end