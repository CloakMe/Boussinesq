%return;
clear;clc;

%btString = '1';
%cString = '90';
%hString = '40';
%orderString = '4';
% additionalInfo == 3 plot and compare Integrals
% additionalInfo == 0 compare solutions
additionalInfo = 3;
if( true )
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '40' ,'2', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '20' ,'2', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '10' ,'2', additionalInfo);
    
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '40' ,'4', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '20' ,'4', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '10' ,'4', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '40' ,'6', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '20' ,'6', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '10' ,'6', additionalInfo);
else
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '20' ,'2', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '10' ,'2', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '05' ,'2', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '20' ,'4', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '10' ,'4', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '05' ,'4', additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '20' ,'6', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '10' ,'6', additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '05' ,'6', additionalInfo);
end