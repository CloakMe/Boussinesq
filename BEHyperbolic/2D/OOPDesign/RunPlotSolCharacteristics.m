clear;
addpath('..\..\..\Common');
btString = { '1' };
cString = { '90' };
hString = { '40' }; %, '10', '20'
orderString = { '6'}; % , '4', '6'
domainLenString = { '512', '256', '128' }; %30 128
solTypeString = { 'Taylor' }; % Taylor  EnergySave
bndString = { 'ZeroBoundary' }; %ZeroBoundary ZeroBnd
%legendString = { 'Taylor O(h^2+{\tau}^2)', 'Taylor O(h^4+{\tau}^4)', 'Taylor O(h^6+{\tau}^6)' }; %'Cons. Scheme O(h^2+{\tau}^2)' 
%legendString = { 'h = 0.05, p=6', 'h = 0.1, p=6', 'h = 0.2, p=6', 'h = 0.05, p=2', 'h = 0.1, p=2', 'h = 0.2, p=2'};
%legendString = { '\Omega_h = [-120,120]x[-108,108]', '\Omega_h = [-60,60]x[-54,54]', '\Omega_h = [-30,30]x[-27,27]'};
legendString = { '\Omega_h = [-512,512]x[-232,232]', '\Omega_h = [-256,256]x[-116,116]', '\Omega_h = [-128,128]x[-58,58]'};
additionalInfo = 'mass';

PlotSolCharacteristics(btString, cString, hString ,orderString, domainLenString, solTypeString, bndString, legendString, additionalInfo);
return;

legendString = { 'Taylor O(h^2+{\tau}^2)', 'Taylor O(h^4+{\tau}^4)', 'Taylor O(h^6+{\tau}^6)',  'Cons. Scheme O(h^2+{\tau}^2)'  }; %
legend(legendString); 


figure(20)
hold on;
plot(tt(1:799:end),EN(1:799:end),'bo' )
hold off;
xlabel('t','FontSize',20);  ylabel('E_h','FontSize',20);
set(gca,'FontSize',16);

figure(20)
hold on;
[indeces, shift] = BEUtilities.GetCommonIndexArray( t, II );
indeces(1) = [];
plot( t(indeces(1:500:end)) ,II(indeces(1:500:end)),'bo' )
hold off;
xlabel('t','FontSize',20);  ylabel('E_h','FontSize',20);
set(gca,'FontSize',18);