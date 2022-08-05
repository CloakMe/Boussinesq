clear;
sx = 10;
%AD = zeros(4,8);
stencil = -3:3;
fdsize = length(GetFiniteDifferenceCoeff( stencil, 2));
fgSizeOrg = fdsize;
AB = zeros(sx);
mid = (fdsize - 1)/2 ;

for i = 1:sx - fdsize + 1
    AB(i+mid,i:i+fdsize-1) = GetFiniteDifferenceCoeff( stencil, 2)';
end

j = 1;
for i=2:-1:-2
    if(i==1)
        c_i = GetFiniteDifferenceCoeff(-fgSizeOrg:0,2)';
        d_i = GetFiniteDifferenceCoeff(-fgSizeOrg+1:1,2)';
        d0_c0 = d_i(1)/c_i(1);
        coeffs = fliplr(d_i(2:end)-d0_c0*c_i(2:end)); 
    elseif(i==0)
        c_i = GetFiniteDifferenceCoeff(-fgSizeOrg:0,2)';
        d_i = GetFiniteDifferenceCoeff(-fgSizeOrg+2:2,2)';
        d0_c0 = d_i(1)/c_i(1);
        coeffs = fliplr(d_i(2:end)-d0_c0*c_i(2:end));   
    else
        stencil = i - 2:i+5;
        coeffs = GetFiniteDifferenceCoeff( stencil, 2)';
    end
    fdsize = length(coeffs);
    AB(j,1:fdsize) = coeffs;
    j = j + 1;
end

j = sx;
for i=2:-1:-2
    if(i==1)
        c_i = GetFiniteDifferenceCoeff(-fgSizeOrg:0,2)';
        d_i = GetFiniteDifferenceCoeff(-fgSizeOrg+1:1,2)';
        d0_c0 = d_i(1)/c_i(1);
        coeffs = fliplr(d_i(2:end)-d0_c0*c_i(2:end)); 
    elseif(i==0)
        c_i = GetFiniteDifferenceCoeff(-fgSizeOrg:0,2)';
        d_i = GetFiniteDifferenceCoeff(-fgSizeOrg+2:2,2)';
        d0_c0 = d_i(1)/c_i(1);
        coeffs = fliplr(d_i(2:end)-d0_c0*c_i(2:end));   
    else
        stencil = i - 2:i+5;
        coeffs = GetFiniteDifferenceCoeff( stencil, 2)';
    end
    fdsize = length(coeffs);
    AB(j,end-fdsize+1:end) = fliplr(coeffs);
    j = j - 1;
end
A = AB(2:end-1,2:end-1)

return;
[ S, D ] = eig( -(Delta_hy'+Delta_hy)/2 );
max(max(diag(D)))
min(min(diag(D)))

Sym_hy = BEUtilities.GetFinDiffMat(this.sy,this.order);
[ SSym, DSym ] = eig( Sym_hy );
max(max(diag(DSym)))
min(min(diag(DSym)))

figure(1)
plot(1:this.sy, diag(DSym), 'g', 1:this.sy, diag(D), 'k', 1:this.sy, diag(D) + diag(DSym), 'r')

Diff_hy = Sym_hy - (Delta_hy'+Delta_hy)/2 ; %+ 2*diag(ones(this.sy,1))
[ SDif, DDif ] = eig( Diff_hy );
max(max(diag(DSum)))
min(min(diag(DSum)))

figure(2)
plot(1:this.sy, diag(DDif), 'b', 1:this.sy, diag(-DSym), 'g', 1:this.sy, diag(DSum - DSym), 'r');
legend('DDif','-DSym', 'DDif-DSym');
Sum_hy(1:10,1:10)

