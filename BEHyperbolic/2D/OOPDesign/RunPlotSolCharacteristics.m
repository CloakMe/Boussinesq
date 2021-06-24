clear
btString = { '1' };
cString = { '90' };
hString = { '10', '20', '40' };
orderString = { '2' };
domainLenString = { '40x80' };
solTypeString = { 'Taylor' };
bndString = { 'WithBoundary' };
legendString = { 'O(h^2+{\tau}^2)', 'O(h^4+{\tau}^4)', 'O(h^6+{\tau}^6)' };
additionalInfo = 'maximum';

PlotSolCharacteristics(btString, cString, hString ,orderString, domainLenString, solTypeString, bndString, legendString, additionalInfo);