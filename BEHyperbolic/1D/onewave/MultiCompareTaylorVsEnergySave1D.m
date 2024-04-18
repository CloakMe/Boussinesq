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
useSecondOrderOnly = 1;
for i = 1:50
    figNumber = i;
    if(  ishandle(i) == true )
        close(figure(i));
    end
end

additionalInfo = 3;
%domainLen = '30'; %128, 30

if( true )
    domainLen = '60'; %128, 30
    CompareTaylorVsEnergySave1D('0000250', '40' ,domainLen, additionalInfo);
    CompareTaylorVsEnergySave1D('0000250', '20' ,domainLen, additionalInfo);
    %CompareTaylorVsEnergySave1D('1', '90', '10' ,domainLen);
else
    domainLen = '30'; %128, 30
    CompareTaylorVsEnergySave1D('3', '45', '20' ,'2', domainLen, additionalInfo);
    CompareTaylorVsEnergySave1D('3', '45', '10' ,'2', domainLen, additionalInfo);
    CompareTaylorVsEnergySave1D('3', '45', '05' ,'2', domainLen, additionalInfo);
    if(useSecondOrderOnly == 1)
        return;
    end
    fprintf('=========================\n'); 
    fprintf('=========================\n\n');
    CompareTaylorVsEnergySave1D('3', '45', '20' ,'4', domainLen, additionalInfo);
    CompareTaylorVsEnergySave1D('3', '45', '10' ,'4', domainLen, additionalInfo);
    CompareTaylorVsEnergySave1D('3', '45', '05' ,'4', domainLen, additionalInfo);

end