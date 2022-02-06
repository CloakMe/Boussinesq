function [maxvalv, max_err, min_err] = xtrct_prop(x,tt,v1,vl)
%maxvalv is the maximum value of the solution along the last time layer
%
%max_err is the maximum of the difference between exact and approximate
%solution along the last time layer
%
%min_err is the minimum of the difference between exact and approximate
%solution along the last time layer

maxvalv = 1.05*max(v1);
    if(isnan(maxvalv))
        maxvalv = 1.5*max(v1);
        if(isnan(maxvalv))
            maxvalv = 5;
        end
    end

    Tr = u_ex(x,tt(end),2);% + u_ex(x-5,t(j),-1.5);
    max_err = max(Tr' - vl);
    min_err = min(Tr' - vl);
    if(max_err == 0 && min_err==0)
        max_err = 0.00001; min_err = -0.00001; 
    end