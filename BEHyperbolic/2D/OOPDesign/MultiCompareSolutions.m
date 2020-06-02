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

for i = 1:50
    figNumber = i;
    if(  ishandle(i) == true )
        close(figure(i));
    end
end
additionalInfo = 4;
if( true )
    CompareTaylorVsEnergySave('1', '90', '40' ,'2', additionalInfo);
    CompareTaylorVsEnergySave('1', '90', '20' ,'2', additionalInfo);
    CompareTaylorVsEnergySave('1', '90', '10' ,'2', additionalInfo);
    
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorVsEnergySave('1', '90', '40' ,'4', additionalInfo);
    CompareTaylorVsEnergySave('1', '90', '20' ,'4', additionalInfo);
    CompareTaylorVsEnergySave('1', '90', '10' ,'4', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorVsEnergySave('1', '90', '40' ,'6', additionalInfo);
    CompareTaylorVsEnergySave('1', '90', '20' ,'6', additionalInfo);
    CompareTaylorVsEnergySave('1', '90', '10' ,'6', additionalInfo);
else
    CompareTaylorVsEnergySave('3', '52', '20' ,'2', additionalInfo);
    CompareTaylorVsEnergySave('3', '52', '10' ,'2', additionalInfo);
    CompareTaylorVsEnergySave('3', '52', '05' ,'2', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorVsEnergySave('3', '52', '20' ,'4', additionalInfo);
    CompareTaylorVsEnergySave('3', '52', '10' ,'4', additionalInfo);
    CompareTaylorVsEnergySave('3', '52', '05' ,'4', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorVsEnergySave('3', '52', '20' ,'6', additionalInfo);
    CompareTaylorVsEnergySave('3', '52', '10' ,'6', additionalInfo);
    CompareTaylorVsEnergySave('3', '52', '05' ,'6', additionalInfo);
end