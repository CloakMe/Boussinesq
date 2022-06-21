function [ name ] = GetICName( ICSwitch )
%UNTITLED Summary of this function goes here
%   Type of Initial Data 
%   ICSwitch == 1 --> Natali Inital approximation (ground state)
%   ICSwitch == 0 --> Christov Inital approximation (
    if(ICSwitch == 1)
       name = 'Natali'; 
       return;
    end
    
    if(ICSwitch ~= 1)
       name = 'Christov'; 
       return;
    end

end

