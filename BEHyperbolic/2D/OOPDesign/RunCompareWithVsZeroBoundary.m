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

additionalInfo = 2; % 
CompareWithVsZeroBoundary('3', '45', '50' ,'4', '20', additionalInfo ) 
CompareWithVsZeroBoundary('3', '45', '50' ,'4', '40', additionalInfo ) 
CompareWithVsZeroBoundary('3', '45', '50' ,'4', '80', additionalInfo ) 
CompareWithVsZeroBoundary('3', '45', '50' ,'4', '160', additionalInfo ) 

return;
legendString = { 'h=0.4, p=4', 'h=0.2, p=4', 'h=0.1, p=4', 'h=0.4, p=6', 'h=0.2, p=6', 'h=0.1, p=6'}; %
legend(legendString); 

