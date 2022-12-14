%return;
clear;clc;

%btString = '1';
%cString = '90';
%hString = '40';
%orderString = '4';
% additionalInfo == 1 compare solutions
% additionalInfo == 2 compare solutions with additional graphs
% additionalInfo == 3 compare integrals
% additionalInfo == 4 compare integrals with additional graphs
% useSecondOrderOnly = 1;
% for i = 1:50
%     figNumber = i;
%     if(  ishandle(i) == true )
%         close(figure(i));
%     end
% end

additionalInfo = 1;
domainLen = '50'; %128, 30, 50
% legendString = { 'h=0.4, {\tau} = 0.04', 'h=0.2, {\tau} = 0.02', 'h=0.1, {\tau} = 0.01'};
% legend(legendString);
tFix = false;
if( strcmp(domainLen,'50') == 1 )
    CompareStartEndSolitons('3', '30', '40' ,'4', domainLen, additionalInfo, 'ro', tFix ) 
    CompareStartEndSolitons('3', '30', '20' ,'4', domainLen, additionalInfo, 'go' ) 
    CompareStartEndSolitons('3', '30', '10' ,'4', domainLen, additionalInfo, 'bo' ) 

    CompareStartEndSolitons('3', '30', '40' ,'6', domainLen, additionalInfo, 'r', tFix ) 
    CompareStartEndSolitons('3', '30', '20' ,'6', domainLen, additionalInfo, 'g' ) 
    CompareStartEndSolitons('3', '30', '10' ,'6', domainLen, additionalInfo, 'b' ) 
    legendString = { 'h=0.4, p=4', 'h=0.2, p=4', 'h=0.1, p=4', 'h=0.4, p=6', 'h=0.2, p=6', 'h=0.1, p=6'}; %
    legend(legendString); 
end

if( strcmp(domainLen,'128') == 1 )
    CompareStartEndSolitons('1', '90', '40' ,'2', domainLen, additionalInfo, 'b', tFix ) 
    CompareStartEndSolitons('1', '90', '20' ,'2', domainLen, additionalInfo, 'gx' ) 
    CompareStartEndSolitons('1', '90', '10' ,'2', domainLen, additionalInfo, 'md' ) 

    CompareStartEndSolitons('1', '90', '40' ,'4', domainLen, additionalInfo, 'go', tFix ) 
    CompareStartEndSolitons('1', '90', '20' ,'4', domainLen, additionalInfo, 'gx' ) 
    CompareStartEndSolitons('1', '90', '10' ,'4', domainLen, additionalInfo, 'gd' ) 
    
    CompareStartEndSolitons('1', '90', '40' ,'6', domainLen, additionalInfo, 'mo', tFix ) 
    CompareStartEndSolitons('1', '90', '20' ,'6', domainLen, additionalInfo, 'mx' ) 
    CompareStartEndSolitons('1', '90', '10' ,'6', domainLen, additionalInfo, 'md' ) 
end

if( strcmp(domainLen,'30') == 1 )
    CompareStartEndSolitons('3', '45', '20' ,'2', domainLen, additionalInfo, 'ko', tFix ) 
    CompareStartEndSolitons('3', '45', '10' ,'2', domainLen, additionalInfo, 'kx' ) 
    CompareStartEndSolitons('3', '45', '05' ,'2', domainLen, additionalInfo, 'kd' ) 
    
    CompareStartEndSolitons('3', '45', '20' ,'4', domainLen, additionalInfo, 'go', tFix ) 
    CompareStartEndSolitons('3', '45', '10' ,'4', domainLen, additionalInfo, 'gx' ) 
    CompareStartEndSolitons('3', '45', '05' ,'4', domainLen, additionalInfo, 'gd' ) 
    
    CompareStartEndSolitons('3', '45', '20' ,'6', domainLen, additionalInfo, 'mo', tFix ) 
    CompareStartEndSolitons('3', '45', '10' ,'6', domainLen, additionalInfo, 'mx' ) 
    CompareStartEndSolitons('3', '45', '05' ,'6', domainLen, additionalInfo, 'md' )     
end

