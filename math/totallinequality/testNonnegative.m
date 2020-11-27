% %°¸Àı1
%  e=0.0001;
%  m1=500;
% m2=500;
% n=700;
% A1=sprand(m1,n,0.1,1/100);
% A2=sprand(m2,n,0.1,1/100);
% b1=rand(m1,1);
% b2=rand(m2,1);
% A=[A1;-A2];
% b=[b1;-b2];
% x0=ones(n,1);
% det=ones(n,1);
% A2(A2<0)=0.5;
% A1(A1<0)=1;
% %°¸Àı2
% A=[1 2 3;4 5 6;-7 -8 -9];
% b=[10 11 -12]';
% x0=[1 1 1]';
% e=0.00001;
% det=[1 1 1]';
% 
% %
addpath('./IPG')
addpath('./exact')
clear 
clc
m1 = 1000; 
m2 = 1000; n = 300;
A1=sprand(m1,n,0.1,1/100);
A2=sprand(m2,n,0.1,1/100);
b1=rand(m1,1);
b2=rand(m2,1);
A=[A1;-A2];
b=[b1;-b2];
x0=ones(n,1);
det=ones(n,1);
[x0, resvec, arvec, tf]=fixedMatrix(A,b,x0,20,1e-15);
tf
% A2(A2<0)=0.5;
% A1(A1<0)=1;
% A=[2,1;2,-1];
% b=[5;-3];
% [x,fk]=fsearchx(A,b,[1;1],0.000001,0.0001)
