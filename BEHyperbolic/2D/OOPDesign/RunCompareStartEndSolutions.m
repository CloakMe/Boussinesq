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

additionalInfo = 2;

if( false )
    domainLen = '128'; %128, 30
    CompareStartEndSolitons('1', '90', '40' ,'6', domainLen, additionalInfo, 'b' ) 
    CompareStartEndSolitons('1', '90', '20' ,'6', domainLen, additionalInfo, 'gx' ) 
    CompareStartEndSolitons('1', '90', '10' ,'6', domainLen, additionalInfo, 'md' ) 
    if( additionalInfo == 2 )
        %legendString = { 'h=0.4, {\tau} = 0.002', 'h=0.2, {\tau} = 0.001', 'h=0.1, {\tau} = 0.0005'}; %
        %legendString = { 'h=0.4, {\tau} = 0.04', 'h=0.2, {\tau} = 0.02', 'h=0.1, {\tau} = 0.01'}; %
        legendString = { 'h=0.4, {\tau} = 0.04', 'h=0.2, {\tau} = 0.02', 'h=0.1, {\tau} = 0.01'}; %
        legend(legendString); 
    end
    return;
    CompareStartEndSolitons('1', '90', '40' ,'4', domainLen, additionalInfo, 'go' ) 
    CompareStartEndSolitons('1', '90', '20' ,'4', domainLen, additionalInfo, 'gx' ) 
    CompareStartEndSolitons('1', '90', '10' ,'4', domainLen, additionalInfo, 'gd' ) 
    
    CompareStartEndSolitons('1', '90', '40' ,'6', domainLen, additionalInfo, 'mo' ) 
    CompareStartEndSolitons('1', '90', '20' ,'6', domainLen, additionalInfo, 'mx' ) 
    CompareStartEndSolitons('1', '90', '10' ,'6', domainLen, additionalInfo, 'md' ) 
    

else
    domainLen = '30'; %128, 30
    CompareStartEndSolitons('3', '45', '20' ,'2', domainLen, additionalInfo, 'ko' ) 
    CompareStartEndSolitons('3', '45', '10' ,'2', domainLen, additionalInfo, 'kx' ) 
    CompareStartEndSolitons('3', '45', '05' ,'2', domainLen, additionalInfo, 'kd' ) 
    
    CompareStartEndSolitons('3', '45', '20' ,'4', domainLen, additionalInfo, 'go' ) 
    CompareStartEndSolitons('3', '45', '10' ,'4', domainLen, additionalInfo, 'gx' ) 
    CompareStartEndSolitons('3', '45', '05' ,'4', domainLen, additionalInfo, 'gd' ) 
    
    CompareStartEndSolitons('3', '45', '20' ,'6', domainLen, additionalInfo, 'mo' ) 
    CompareStartEndSolitons('3', '45', '10' ,'6', domainLen, additionalInfo, 'mx' ) 
    CompareStartEndSolitons('3', '45', '05' ,'6', domainLen, additionalInfo, 'md' )     
end

