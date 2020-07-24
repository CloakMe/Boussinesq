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
useSecondOrderOnly = true;

for i = 1:10
    figNumber = i;
    if(  ishandle(i) == true )
        close(figure(i));
    end
end
additionalInfo = 2;
domainLen = '80';
if( true )
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '40' ,'2', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '20' ,'2', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '10' ,'2', domainLen, additionalInfo);
    if(useSecondOrderOnly == 1)
        return;
    end
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '40' ,'4', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '20' ,'4', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '10' ,'4', domainLen, additionalInfo);

    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '40' ,'6', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '20' ,'6', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('1', '90', '10' ,'6', domainLen, additionalInfo);
else
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '20' ,'2', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '10' ,'2', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '05' ,'2', domainLen, additionalInfo);
    if(useSecondOrderOnly == 1)
        return;
    end
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '20' ,'4', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '10' ,'4', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '05' ,'4', domainLen, additionalInfo);
    
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '20' ,'6', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '10' ,'6', domainLen, additionalInfo);
    CompareTaylorWithBoundaryVsZeroBoundary('3', '52', '05' ,'6', domainLen, additionalInfo);
end