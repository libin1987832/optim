% compare simple iterator and newton method 
A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];

[q,r]=qr(A);

x=[-1;0.05;0];
y=[1;0.05:0];


A=[1,1;-1,-1;1,0;6,3];
b=[1;1;-0.5;-2];
% output
d=lineData(A,b,[-2,1],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
line([-2/3,1],[2/3,-1],'LineStyle','--');


