clear;
btString = '1';
cString = '90';
hString = '40';
orderString = '4';
[x1,y1,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString );
[x2,y2,uEnTaylor] = GetBEEngineTaylorSol( btString, cString, hString, orderString );

if( length( x1 ) ~= length( x2 ) || length( y1 ) ~= length( y2 ) )
    fprintf('Different sizes in X, Y - dimensions!');
elseif( sum( x1 ~= x2 ) && sum( y1 ~= y2 ) )
    fprintf('Different values in x, y dimension vectors!');
else
    fprintf('sx = %d, sy = %d', length(x1), length(y1));
end

viewTypeX = 0;
viewTypeY = 90;
  
figure(14)
mesh(x1,y1,uEnSave');
view( viewTypeX, viewTypeY );
colorbar;
title('solution Energy Save');
xlabel('x');            ylabel('y');

figure(15)
mesh(x2,y2,uEnTaylor');
view( viewTypeX, viewTypeY );
colorbar;
title('solution Taylor');
xlabel('x');            ylabel('y');

figure(16)
mesh(x1,y1,(uEnTaylor- uEnSave)');
view( viewTypeX, viewTypeY );
colorbar;
title('solution difference');
xlabel('x');            ylabel('y');

