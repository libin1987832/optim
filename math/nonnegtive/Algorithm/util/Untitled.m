a=-1:0.01:2;
A=[1,2;3,4];
b=[1;1];
d=[2;3];
x0=[0;0];
y=[];
for i=a
    f=dFM(A,b,x0+i*d);
    y=[y;f'*f];
end
plot(a,y);