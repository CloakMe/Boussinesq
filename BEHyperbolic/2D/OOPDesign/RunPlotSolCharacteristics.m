clear
btString = { '3' };
cString = { '45' };
hString = { '05', '10', '20' };
orderString = { '2' };
domainLenString = { '30' };
solTypeString = { 'EnergySave' };
bndString = { 'ZeroBoundary' };
%legendString = { 'Taylor O(h^2+{\tau}^2)', 'Taylor O(h^4+{\tau}^4)', 'Taylor O(h^6+{\tau}^6)' }; Conservative Scheme O(h^2+{\tau}^2)
legendString = { 'h = 0.05', 'h = 0.1', 'h = 0.2'};
additionalInfo = 'energy';

PlotSolCharacteristics(btString, cString, hString ,orderString, domainLenString, solTypeString, bndString, legendString, additionalInfo);
return;

legendString = { 'Taylor O(h^2+{\tau}^2)', 'Taylor O(h^4+{\tau}^4)', 'Taylor O(h^6+{\tau}^6)',  'Conservative Scheme O(h^2+{\tau}^2)'  }; %
legend(legendString); 


figure(20)
hold on;
plot(tt(1:799:end),EN(1:799:end),'bo' )
hold off;
title('Energy functional');
xlabel('time "t"');  ylabel('EN');

figure(20)
hold on;
[indeces, shift] = BEUtilities.GetCommonIndexArray( t, II );
indeces(1) = [];
plot( t(indeces(1:399:end)) ,II(indeces(1:399:end)),'bo' )
hold off;
title('Integral');
xlabel('time "t"');  ylabel('Integral');