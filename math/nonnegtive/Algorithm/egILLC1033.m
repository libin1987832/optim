addpath('FM');
addpath('util');
[A,rows,cols,entries,rep,field,symm]=mmread('./util/illc1033.mtx');
b=ones(rows,1);
% x0=zeros(cols,1)+1.23;
% x0=zeros(cols,1)+1.13; 结果差了
% x0=zeros(cols,1)+1.5; 我们结果好 很明显
% [m,n]=size(A);
[xk1,fk1,xkArr1]=hybrid1(x0,A,b);
[xk2,fk2,xkArr2]=hybrid2(x0,A,b);