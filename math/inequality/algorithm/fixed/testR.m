% ����ϵͳ������ R\b ���Լ���ѭ��д �����Լ�д�ı�ϵͳ����Ҫ�� 
% m= 100 b=50 0.000169 our 0.0018(first)  0.000097 0.00004(many test)
addpath('../../dataInequality')
addpath('../../algorithmInequality')
clear
clc
m = 10000;
n = 5000;
[A,b,x00] = randInequality(m,n,2,-1);
%tic
[Q,R] = qr(A);
bx=b;
tic
xr=R\bx;
toc
tic
x=x00;
for i=n:-1:1
    for j=i:n
        bx(i)=bx(i)-R(i,j)*x(j);
    end
    x(i)=bx(i)/R(i,i);
end
toc