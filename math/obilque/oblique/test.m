clc;clear;close all;
%初始值
m=3;n=5;
A=CreatA(m,n);
b=rand(m,1);RRE=0.5e-6;x_0=zeros(n,1);

if m>=n
    x_exact=(A'*A)\A'*b;
    b_A=b-A*x_exact;
else
    b_A=b-A'\(A'*b);
end
%%
[x,k,R]=GSO_improve(A,b,n,RRE,x_0,b_A);
disp('GSO迭代次数:');
disp(k);
disp(R);
subplot(1,3,1);
plot3(R(1,:),R(2,:),R(3,:))
subplot(1,3,2);
[x,y]=meshgrid(-2:0.1:1);
for j=1:1:n
    z=(A(1,j)*x+A(2,j)*y)/A(3,j);
    mesh(x,y,z);
    hold on;
end



