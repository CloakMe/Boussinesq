% constants
x_st =  -15.0;   y_st = -15.0;
x_end = 15.0;   y_end = 22.0;
x_st2 = -1.0;   y_st2 = -1.0;
x_end2 = -1.0;   y_end2 = -1.0;
compBox = struct('x_st',{x_st},'x_end',{x_end},'y_st',{y_st},'y_end',...
    {y_end},'x_st2',{x_st2},'x_end2',{x_end2},'y_st2',{y_st2},'y_end2',{y_end2});

h = 0.05; 
x=x_st:h:x_end; y=y_st:h:y_end;
al = -1;%99979 izb
beta1 = 3;
beta2 = 1;
c = 0.55;

vc=1;
order=2; %O( h^ord + tau^ord )

sy = length(y)
sx = length(x)

ox1 = ones(1,sy);
X = x'*ox1;
if(sx~=sy)
    oy1 = ones(1,sx);
    Y = (y'*oy1)';
else
    Y = (y'*ox1)';
end

