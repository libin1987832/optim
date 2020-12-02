
y1 = @(x,k) 4*sin(pi*x*k)/(pi*k);
x = -1:0.01:1;
for i = 1:10
    yy1=y1(x,1);
end
array
