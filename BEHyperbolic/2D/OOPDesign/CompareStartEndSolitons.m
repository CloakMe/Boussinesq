%c = 0.9
%P:\PhDWork\OneSoliton_bt1\ChristovIC_40_80_bt1_c090_h010_O(h^6).mat
%P:\PhDWork\OneSoliton_bt1\Hyperb_40_80_bt1_c090_h010_O(h^6).mat
%clear
%c = 0.52;
%load "C:\Boussinesq\BEHyperbolic\2D\OOPDesign\SavedWorkspaces\Hyperb_40_bt3_c052_h005_O(h^6).mat"
%load "C:\Boussinesq\BEEliptic\Boussinesq2D\WithBoundary\ChristovIC_40_bt3_c052\Oh6\ChristovIC_40_bt3_c052_h005_O(h^6).mat"

viewTypeX = 90;
viewTypeY = 90;

x_idx = GetIdx( x, 0 )
y_start_idx = GetIdx( y, 0 )
y_end_idx = GetIdx( y, t(end)*c )

y_half_size = min( length(y) - y_end_idx - 20, floor( (y_end_idx - y_start_idx)/2 ) )
x_half_size = length(x) - x_idx - 20

indeces_x = x_idx - x_half_size:x_idx + x_half_size;
start_indeces_y = y_start_idx - y_half_size:y_start_idx + y_half_size;
end_indeces_y = y_end_idx - y_half_size:y_end_idx + y_half_size;

%return;
x_vec = x(indeces_x);
y_start = y(start_indeces_y);
y_end = y(end_indeces_y);

vl_start = bigU( indeces_x, start_indeces_y );
vl_end = vl( indeces_x, end_indeces_y );

figure(1)
mesh(x_vec, y_start, vl_start');
title(['Solution at time: ', num2str(0)]);
xlabel('x'); ylabel('y');
colorbar;
axis tight;
view( viewTypeX, viewTypeY );

figure(2)
mesh(x_vec, y_end, vl_end');
title(['Solution at time: ', num2str(t(end))]);
xlabel('x'); ylabel('y');
colorbar;
axis tight;
view( viewTypeX, viewTypeY );

figure(3)
mesh(x_vec, y_start, (vl_start-vl_end)');
title('Solution difference between start and end time');
xlabel('x'); ylabel('y');
colorbar;
axis tight;
view( viewTypeX, viewTypeY );
return;
figure(4)
mesh(x, y, vl');
title('Solution');
xlabel('x'); ylabel('y');
colorbar;
axis tight;