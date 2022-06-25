%%�����AΪ�Գƣ�����������BΪ��λ���������ֵ��������.
%% �������n�׳�ʼ����A��B
clc
clear
addpath('../')
n=50;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
a=-21;
b=10;
r=a+b.*unidrnd(100,n,1);% �����䣨a��b���ڲ������ȷֲ���nά����
Q=diag(r);
A=C'*Q*C; %��ʼ����A
B=eye(n);%��ʼ����B
%% �����ʼ����
r=zeros(n,1);
for i=1:n
 r(i)=min(A(:,i)*B(i,i)-A(i,i)*B(:,i));
 if r(i)>=0
  disp('�þ�����Ҫ')
 end
end
s=find(r==max(r));
s=s(1);
m=zeros(n,1);
m(s)=1;
x1=m; %��ʼ����x1
%% �������к���


[xk,i,h,lamdab]=SPL(A,B,x1,n);
[iD]=SSQPD(A,x1,n); %����SSQP(D)���������ʱ��͵�������
R=max(abs(eig(A)));    %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(n,1);
e(1)=R+0.01;
I=eye(n);
% M=A;%A�Գ�����
M=A+(R+0.01)*I;
M = diag(diag(M));
sigma0 = 0;
tic;[x,  iter, error] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);toc
lambda = (x' * A * x) / (x' * B * x);
disp(['sqp:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda)) ',iter=' num2str(iter)])
[iI]=SSQPI(A,x1,n); %����SSQP(I)���������ʱ��͵�������
[iG]=SQPG(A,x1,n); %����SQP(G)���������ʱ��͵�������
[i]=bas1(A,x1,n); %����BAS���������ʱ��͵�������
