%% SPL（别人的算法）
function [xk,i,h,lamdab]=SPL(A,B,epsilon,epsilon2,e,n,e0)
 x0=zeros(n,1);
 x0(1)=1;
 tic %计时
 i=0;
 options = optimoptions('quadprog','Display','off');%为了不输出 quadprog中间的过程
 while 1
lamdab=(x0'*A*x0)/(x0'*B*x0);
yk=-(lamdab*B*x0) ;
xk=BBP(A,yk,e,n,epsilon2);
% xk=quadprog(A,yk,[],[],e,1,e0,[],[],options); %调用内置函数quadprog求解xk.
dk=xk-x0;
w=(A-lamdab*B)*xk;
gc=(xk'*w)/2;
if abs(gc)<epsilon
break
end     
x0=xk; %迭代更新
i=i+1;
 end 
 toc %计时结束
w=(A-lamdab*B)*xk;
h1=min(xk,w);
h=norm(h1); % 互补性
end