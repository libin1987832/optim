%% SPL�����˵��㷨��
function [xk,i,h,lamdab]=SPL(A,B,x1,n)
epsilon=10^(-5);          %����
epsilon2=10^(-6);        %������ľ���
e=ones(1,n);%ȫ��1�ĺ�����
x0=x1;
tic %��ʱ
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
x0=xk; %��������
i=i+1;
 XZI5=[XZI5   i];
 end 
 toc %��ʱ����
dlmwrite('XZB5.txt', XZI5, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZB5.txt', YZI5, 'precision', '%5f', 'delimiter', '\t');
lamdab=lamdab
i
A=M;
w=(A-lamdab*B)*xk;
min(w);
h1=min(xk,w);
h=norm(h1); % ������
end