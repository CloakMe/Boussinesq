function [d2vz, d3vz, d4vz, d5vz] = calc_der(vz,dvz,sdh,sIdh,sdh11,sIdh11,h,alpha,beta)
 %==================
    s_isdh = size(sIdh,1)*size(sIdh,2);
    b = alpha*beta*vz.*vz + (beta-1)*vz;
    
    deltab = BMM(sdh,b',sdh11);    deltav  =  BMM(sdh,vz',sdh11);
    if(s_isdh == 3)  y = diag3solv(sIdh,deltab);  
    else  y = pentsolv(sIdh11,sIdh,deltab); 
    end
    d2vz = (y + deltav'/h^2);
    
    %==================    
    b = 2*alpha*beta*dvz.*vz + (beta-1)*dvz;
    
    deltab = BMM(sdh,b',sdh11);    deltav  =  BMM(sdh,dvz',sdh11);
    
    if(s_isdh == 3)  y = diag3solv(sIdh,deltab);  
    else  y = pentsolv(sIdh11,sIdh,deltab); 
    end
    d3vz = (y + deltav'/h^2);
    
    %==================    
    b = 2*alpha*beta*(dvz.*dvz + vz.*d2vz) + (beta-1)*d2vz;
    
    deltab = BMM(sdh,b',sdh11);    deltav  =  BMM(sdh,d2vz',sdh11);
    
    if(s_isdh == 3)  y = diag3solv(sIdh,deltab);  
    else  y = pentsolv(sIdh11,sIdh,deltab); 
    end
    d4vz = (y + deltav'/h^2);
    
    %==================  
    b = 2*alpha*beta*(3*dvz.*d2vz + vz.*d3vz) + (beta-1)*d3vz;

    deltab = BMM(sdh,b',sdh(3)+1);    deltav  =  BMM(sdh,d3vz',sdh11);
    
    if(s_isdh == 3)  y = diag3solv(sIdh,deltab);
    else  y = pentsolv(sIdh11,sIdh,deltab);
    end
    d5vz = (y + deltav'/h^2);

    %==================  