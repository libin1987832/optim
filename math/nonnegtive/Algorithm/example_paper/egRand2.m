clc
clear
addpath('../FM')
% m=20;
% ratio=0.3;
% n=ceil(ratio*m);
%  A=2*rand(m,n)-1;
%  b=2*rand(m,1)-1;
%  save('ddf1.mat','A','b','m','n')
 load('ddf1.mat')

% for m=10:20:50
%     for ratio=0.1:0.3:0.5
% m=300;
% ratio=0.5;
% n=ceil(ratio*m);
% A=2*rand(m,n)-1;
% b=2*rand(m,1)-1;
x0=zeros(n,1);
mIter=300;
[xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b,mIter);
[xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b,mIter);
[xk6,fk1,xkArr1,countF1,countN1]=hybrid6(x0,A,b,10,mIter);
% [[1:m]' b-A*xk1 b-A*xk2 b-A*xk6]'
[xk7,fk1,xkArr1,countF1,countN1]=hybrid7(x0,A,b,1,mIter);

%     end
% end

