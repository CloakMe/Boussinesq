function [d2vz, d3vz, d4vz, d5vz] = calc_der_sub(vz,dvz,sdh,sIdh,sdh11,sIdh11,h,alpha,beta)
 %==================
    intermediate =  beta * vz - 1/h^2 * BMM(sdh,vz',sdh11)' - alpha*beta*vz.*vz;
    d2vz = 1/h^2 * BMM(sdh,intermediate',sdh11)';
    
    %==================    
    intermediate =  beta * vz - 1/h^2 * BMM(sdh,vz',sdh11)' - 2*alpha*beta*dvz.*vz;
    d3vz = 1/h^2 * BMM(sdh,intermediate',sdh11)';
    
    %==================    
    intermediate =  beta * vz - 1/h^2 * BMM(sdh,vz',sdh11)' - 2*alpha*beta*(dvz.*dvz + vz.*d2vz);
    d4vz = 1/h^2 * BMM(sdh,intermediate',sdh11)';
    
    %==================  
    intermediate =  beta * vz - 1/h^2 * BMM(sdh,vz',sdh11)' - 2*alpha*beta*(3*dvz.*d2vz + vz.*d3vz) + (beta-1)*d3vz;
    d5vz = 1/h^2 * BMM(sdh,intermediate',sdh11)';    
    %==================  