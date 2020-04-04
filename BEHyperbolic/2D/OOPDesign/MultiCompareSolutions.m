clear;clc;
%return;
%btString = '1';
%cString = '90';
%hString = '40';
%orderString = '4';
% additionalInfo == 3 plot and compare Integrals
% additionalInfo == 0 compare solutions
additionalInfo = 0;
if( true )
    CompareSolutions('1', '90', '40' ,'2', additionalInfo);
    CompareSolutions('1', '90', '20' ,'2', additionalInfo);
    CompareSolutions('1', '90', '10' ,'2', additionalInfo);
    
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareSolutions('1', '90', '40' ,'4', additionalInfo);
    CompareSolutions('1', '90', '20' ,'4', additionalInfo);
    CompareSolutions('1', '90', '10' ,'4', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareSolutions('1', '90', '40' ,'6', additionalInfo);
    CompareSolutions('1', '90', '20' ,'6', additionalInfo);
    CompareSolutions('1', '90', '10' ,'6', additionalInfo);
else
    CompareSolutions('3', '52', '20' ,'2', additionalInfo);
    CompareSolutions('3', '52', '10' ,'2', additionalInfo);
    CompareSolutions('3', '52', '05' ,'2', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareSolutions('3', '52', '20' ,'4', additionalInfo);
    CompareSolutions('3', '52', '10' ,'4', additionalInfo);
    CompareSolutions('3', '52', '05' ,'4', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareSolutions('3', '52', '20' ,'6', additionalInfo);
    CompareSolutions('3', '52', '10' ,'6', additionalInfo);
    CompareSolutions('3', '52', '05' ,'6', additionalInfo);
end