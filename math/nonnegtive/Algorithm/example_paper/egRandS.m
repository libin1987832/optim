% paper rand [1,1] example
%<<<<<<< HEAD
clc
clear
addpath('../FM')

for m=1000:1000:50000
    for ratio=0.1:0.2:1
        n=ceil(ratio*m);
        A=sprandn(m,n,0.4);
        b=sprandn(m,1,0.4);
        x0=sparse(n,1);
        t=clock;
        [xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b);
        tf1=etime(clock,t);
%         [xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);
t=clock;
        [xk6,fk1,xkArr1,countF6,countN6]=hybrid6(x0,A,b,10);
        tf2=etime(clock,t);
%         [xk7,fk1,xkArr1,countF7,countN7]=hybrid7(x0,A,b,10);
        r=b-A*xk1;
        r(r<0)=0;
        fprintf('$ %d \\times %d $ & %4.2f & %g & %d & %d & %4.2f &',m,n,full(0.5*(r'*r)),full(norm(A'*r)),countF1,countN1,tf1);
        r=b-A*xk6;
        r(r<0)=0;
        fprintf('%4.2f & %g & %d & %d & %4.2f\\\\% 4.2f\n',full(0.5*(r'*r)),full(norm(A'*r)),countF6,countN6,tf2,tf1/tf2);
    end
end
% [xk2,fk2,xkArr2,countF2,countN2]=hybrid4(x0,A,b);

