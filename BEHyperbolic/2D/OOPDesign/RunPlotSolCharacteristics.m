clear
btString = { '3' };
cString = { '45' };
hString = { '05', '10', '20' };
orderString = { '2' };
domainLenString = { '30' };
solTypeString = { 'Taylor' };
bndString = { 'ZeroBoundary' };
legendString = { 'O(h^2+{\tau}^2)', 'O(h^4+{\tau}^4)', 'O(h^6+{\tau}^6)' };
additionalInfo = 'maximum';

PlotSolCharacteristics(btString, cString, hString ,orderString, domainLenString, solTypeString, bndString, legendString, additionalInfo);