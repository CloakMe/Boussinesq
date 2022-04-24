function [res] = LE(vkmo,vk,vkpo,sdh,sIdh,sdh11,sIdh11,h,tau,sgm)
   vt = (vkpo - vkmo)/(tau*2);
   sumv = (vk+ vkpo);
   svt = size(vt,1);
            
   Idhvt = BMM(sIdh,vt',sIdh11);     dhsumv = BMM(-sdh,sumv',-sdh11);
   
   s_isdh = size(sIdh,1)*size(sIdh,2);
   if(s_isdh == 3)  vec1 = diag3solv(-sdh/h^2,vt);
   else  vec1 = pentsolv(-sdh11/h^2,-sdh/h^2,vt);
   end
       res=h*(vec1')*vt + h*(vt')*vt + tau^2*(sgm-1/4)*(Idhvt*vt)/h +...  
           h*(( sumv + dhsumv'/h^2)')*sumv/4;
       