function plotFM(A,b,x0,x1)
for x=x0:0.01:x1
    y=b-A*x;
    y(y<0)=0;
    ny=0.5*y'*y;
    plot(x,ny,'.');
    hold on
end
