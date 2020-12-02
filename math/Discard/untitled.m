clear 
clc
y1 = @(x,k) 2*(sin(pi*k)/(pi*k)-cos(pi*k))*sin(pi*x*k)/(pi*k);
x = -3:0.01:3;
ss= [];
for i = 15:20
    yy1 = zeros(size(x));
    for j = 1:20
        if j <i
            yy1=yy1+y1(x,j);
        end    
    end
    ss(i).Y = yy1;
    ss(i).X = x;
end
figure 
hold on
arrayfun(@(a) plot(a.X,a.Y),ss)
hold off