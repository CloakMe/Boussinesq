clear;

aprErr(1,:) = '2';
aprErr(2,:) = '4';
aprErr(3,:) = '6';
aprErr = cellstr(aprErr);

aa = 20;
if (aa > 30)
    stepSize(1,:) = '40';
    stepSize(2,:) = '20';
    stepSize(3,:) = '10';
else
    stepSize(1,:) = '20';
    stepSize(2,:) = '10';
    stepSize(3,:) = '05';
end
stepSize = cellstr(stepSize);

xVec = 2:2:6;
yVec = xVec;
for il = 1:3
    for jl = 1:3        
        name = strcat('SavedWorkspaces\Hyperb_40_bt3_c052_h0', stepSize(il), '_O(h^', aprErr(jl), ').mat' );
        warning('off','all');
        load ( name{1} );
        warning('on','all');
        ggg =6;
        yVec(jl) = elapsedTime/60;
        %clear;
    end
    stepSize(il)
    yVec = yVec
end

