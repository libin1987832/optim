%%�����AΪ�Գƣ�����������BΪ�Գ��������������ֵ��������.
%% �������n�׳�ʼ����A��B
clc
clear
n=100;
C=importdata('C2.txt');
a=100;
b=1;
z1=0:(a-b)/(n-1):1;
z1=z1+b;
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
Q=diag(z);% �����䣨a��b���ڲ������ȷֲ���nά����
A=C'*Q*C; %��ʼ����A

I=eye(n);
B=importdata('B3.txt'); %��ʼ����B
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
% x1=unidrnd(2,n,1)
%% �������к���
% [iG]=SQPGB(A,x1,B,n,I); %����SQP(G)���������ʱ��͵�������
% M=A;%A�Գ�����
% O=inv(B)*A;
% R=max(abs(eig(O))); %A���װ뾶  
% sigma0 = R + 0.01;
% [x, error] = sqpEicp(A, B, M, x1, sigma0, 1e-5, 1e-5, 10000, 0);
[iD]=SSQPDB(A,x1,B,n,I); %����SSQP(D)���������ʱ��͵�������
M=diag(A);%A�Գ�����
O=inv(B)*A;
R=max(abs(eig(O))); %A���װ뾶  
sigma0 = R + 0.01;
R1=max(abs(eig(A)))+0.01;
 M=A+R1*I;%A�ԳƷ�����
%  M=A; %A�Գ�����
M=diag(diag(M));
[x, error] = sqpEicp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-12, 10000, 0);
% [i]=bas1B(A,x1,B,n); %����BAS���������ʱ��͵�������
% [iI]=SSQPIB(A,x1,B,n); %����SSQP(I)���������ʱ��͵�������
% [xk,i,h,lamdab]=SPL(A,B,x1,n);