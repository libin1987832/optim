%%�����AΪ�Գƣ�����������BΪ�Գ��������������ֵ��������.
%% �������n�׳�ʼ����A��B
clc
clear
n=10;
C=unidrnd(10,n,n);
%C = load('C');
%[C,Y]=qr(C.C,0);
[C,Y]=qr(C,0);
a=100;
b=1;
z1=0:(a-b)/(n-1):1;
z1=z1+b;
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
I=eye(n);
Q=diag(z);% �����䣨a��b���ڲ������ȷֲ���nά����
A=C'*Q*C; %��ʼ����A
B=I+C'*C; %importdata('B3.txt'); %��ʼ����B

%  n = 4;
 n = 6;
I=eye(n);
%  A=[5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
A = [0 1 0 1 0 0;1 0 0 1 0 0;0 0 0 1 0 0;1 1 1 0 1 1;0 0 0 1 0 1;0 0 0 1 1 0] ;
 B = I;
M = I*3;
%  M=A; %A�Գ�����
epsx = 0;
epsxlambda = -1e-3;
x1 = ones(n,1)/n;
O=inv(B)*A;
R=max(abs(eig(O))); %A���װ뾶  
sigma0 = R + 0.01;
[xa, iter, error] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-8, 1000, 0);
nx = size(xa,2);
xd = []
for i = 1 : nx
    x = xa(:,i) ./ sum(xa(:,i));
lambda = (x' * A * x) / (x' * B * x);
%xd = [xd [lambda; x]];
xd = [xd lambda];
%disp(['FQP_QUAD:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda))])
end
 sort(uniquetol(xd,1e-3),'descend')'