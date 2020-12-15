clear 
clc
n = 3;
M = [1,0,0;2,1,0;2,2,1];
q = [-8;-12;-14];
x0 = zeros(n,1) + 1;
[w,z,retcode] = LCPSolveL(M,q)

% M = [10,0,-2;2,0.1,-0.4;0,0.2,0.1];
% q = [10;1;-1];
% Mz + q
[w2,z,retcode] = Murty(M,q);
[w3,z,retcode] = Graves(M,q);
[w1,z,retcode] = Bard(M,q);

xk=MPRGP(M,q,x0,0.1,0.01,0.02,zeros(n,1),30);
xk
[w';w1';w2';w3']


