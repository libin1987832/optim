                                                                                                                                                                                     %%�����AΪ�Գƣ�����������BΪ�Գ��������������ֵ��������.
%% �������n�׳�ʼ����A��B
clc
clear
n=200;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
a=1;
b=1;
z1=a+b.*unidrnd(100,n,1);
Q=diag(z1);
A=C'*Q*C; %��ʼ����A

C=-8+22*rand(n,n);
I=eye(n);
B=C'*C+I;%��ʼ����BC'*C+

%% �����ʼ����
r=zeros(n,1);
for i=1:n
 r(i)=min(A(:,i)*B(i,i)-A(i,i)*B(:,i));
 if r(i)>=0
  disp('�þ�����Ҫ')
 end
end 
s=find(r==max(r));
m=zeros(n,1);
m(s)=1;
x1=m; %��ʼ����x1
%% �������к���

% [i]=bas1B(A,x1,B,n); %����BAS���������ʱ��͵�������
% 
 [iG]=SQPGB(A,x1,B,n,I); %����SQP(G)���������ʱ��͵�������
  [x, error] = sqpEicp(A, B, M, x0, sigma0, eps, bisectEps, maxIT, debug)
% 
% [iD]=SSQPDB(A,x1,B,n,I); %����SSQP(D)���������ʱ��͵�������
% 
% [iI]=SSQPIB(A,x1,B,n); %����SSQP(I)���������ʱ��͵�������
% % % % 
% [xk,i,lamdab]=SPL(A,B,x1,n);

