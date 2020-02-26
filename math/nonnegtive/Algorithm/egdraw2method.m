% compare simple iterator and newton method 
A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];

[q,r]=qr(A);

x=[-1;0.05;0];
y=[1;0.05:0];

line([0,-1],[0,1])
line([0,-1],[0,1])
line([-0.5,-0.5],[0,1])
line([-1,0],[0,2/3])


