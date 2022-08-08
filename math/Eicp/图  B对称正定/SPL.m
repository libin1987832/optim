%% SPL（别人的算法）
function [xk,i,h,lamdab]=SPL(A,B,x1,n)
epsilon=10^(-5);          %精度
epsilon2=10^(-6);        %子问题的精度
e=ones(1,n);%全是1的横向量
x0=x1;
tic %计时
i=0;
XZI5=[];
YZI5=[];
 M=A;
% A=A+ 25*B;
 while 1
lamdab=(x0'*A*x0)/(x0'*B*x0);
yk=-(lamdab*B*x0) ;
xk=BBP(A,yk,e,n,epsilon2);
dk=xk-x0;
if norm(dk)<=epsilon
break
end  
YZI5=[YZI5   norm(dk)];
x0=xk; %迭代更新
i=i+1;
 XZI5=[XZI5   i];
 end 
 toc %计时结束
dlmwrite('XZB5.txt', XZI5, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZB5.txt', YZI5, 'precision', '%5f', 'delimiter', '\t');
lamdab=lamdab
i
A=M;
w=(A-lamdab*B)*xk;
min(w);
h1=min(xk,w);
h=norm(h1); % 互补性
end