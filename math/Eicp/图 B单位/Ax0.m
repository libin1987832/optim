%%�����AΪ�Գƣ�����������BΪ��λ���������ֵ��������.
%% �������n�׳�ʼ����A��B
clc
clear
n=100;
C=importdata('C1.txt')
z=zeros(n,1);
a=100;
b=1;
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
Q=diag(z);% �����䣨a��b���ڲ������ȷֲ���nά����
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
[iG]=SQPG(A,x1,n); %����SQP(G)���������ʱ��͵�������
[iD]=SSQPD(A,x1,n); %����SSQP(D)���������ʱ��͵�������
[iI]=SSQPI(A,x1,n); %����SSQP(I)���������ʱ��͵�������
[i]=bas1(A,x1,n); %����BAS���������ʱ��͵�������
[xk,i,h,lamdab]=SPL(A,B,x1,n);