C=[1,2,3;2,5,7;3,7,9];
q=[-5;-12;0];
xs=[1;2;0];
x0=[1;1;1];
addpath('../');
addpath('../test');
xka=[];
for i=1:10
[xks,ress]=splitS(C,q,1,x0,1);
x0=xks;
xka=[xka xks];
end
xka