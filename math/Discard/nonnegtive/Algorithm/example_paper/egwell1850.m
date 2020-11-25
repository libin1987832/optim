addpath('FM');
addpath('util');
[A,rows,cols,entries,rep,field,symm]=mmread('./util/well1850.mtx');
b=ones(rows,1);
% x0=zeros(cols,1)+11.1; ²î
 x0=zeros(cols,1)-1045;

% [m,n]=size(A);
[xk1,fk1,xkArr1]=hybrid1(x0,A,b);
[xk2,fk2,xkArr2]=hybrid2(x0,A,b);