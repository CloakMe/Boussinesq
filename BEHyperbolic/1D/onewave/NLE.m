function [res]=NLE(vk,vkp1,sx,h,al,bt)
    %G = inline('(al*bt*vv^3)/3 + ((bt-1)*vv^2)/2','vv','al','bt');
    
%for i=1:sx
    res =  h*((al*bt*(sum(vk.^3)))/3 + (al*bt*(sum(vkp1.^3)))/3)  + h*( ((bt-1)*sum(vk.^2))/2 + ((bt-1)*sum(vkp1.^2))/2);
%end



