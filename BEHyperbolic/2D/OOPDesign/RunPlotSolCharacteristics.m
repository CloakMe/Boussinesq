clear
btString = { '1' };
cString = { '90' };
hString = { '10', '20', '40' };
orderString = { '2' };
domainLenString = { '128' };
solTypeString = { 'Taylor' };
bndString = { 'ZeroBoundary' };
legendString = { 'Taylor O(h^2+{\tau}^2)', 'Taylor O(h^4+{\tau}^4)', 'Taylor O(h^6+{\tau}^6)' }; %,  'Conservative Scheme O(h^2+{\tau}^2)' 

additionalInfo = 'mass';

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