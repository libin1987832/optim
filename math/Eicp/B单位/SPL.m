%% SPL�����˵��㷨��
function [xk,i,h,lamdab]=SPL(A,B,x1,n)
epsilon=10^(-5);          %����
epsilon2=10^(-6);        %������ľ���
e=ones(1,n);%ȫ��1�ĺ�����
x0=x1;
tic %��ʱ
i=0;
M=A;
A=A+ 25*B;
 while 1
lamdab=(x0'*A*x0)/(x0'*B*x0);
yk=-(lamdab*B*x0) ;
xk=BBP(A,yk,e,n,epsilon2);
dk=xk-x0;
if norm(dk)<=epsilon
break
end     
x0=xk; %��������
i=i+1;
 end 
 toc %��ʱ����
 lamdab=lamdab-25
 i
 A=M;
w=(A-lamdab*B)*xk;
min(w);
h1=min(xk,w);
h=norm(h1); % ������
end