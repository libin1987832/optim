addpath('../');
addpath('../../LCP')
clear
clc
[A,b,x0] = readData(2,0,0,0);
[m,n] = size(A);
AI = [A,-speye(m)];
M = AI' * AI;
q = AI'  * b;
[w,z,retcode] = Lemka(M,q,1e-10,20);
% % [s,i,c]=exponent(A,b,1e-10);
% for i = 1:c
% [rpk2, normr2, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, s(:,i) , [], 1);
% fprintf('%d, %g, %g \n',i,normr2,normKKT2);
% end
x = z(1:n);
[rpk2, normr2, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, x , [], 1);