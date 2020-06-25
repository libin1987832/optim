% paper rand [1,1] example
%<<<<<<< HEAD
clc
clear
addpath('../FM')
tfA=[];
dfA=[];
dim=[];
for m=100:100:200
    for ratio=0.1:0.2:1
        n=ceil(ratio*m);
        A=2*rand(m,n)-1;
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
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
        df1=norm(A'*r);
        fprintf('$ %d \\times %d $ & %4.2f & %g & %d & %d & %4.2f &',m,n,0.5*(r'*r),df1,countF1,countN1,tf1);
        r=b-A*xk6;
        r(r<0)=0;
        df6=norm(A'*r);
        fprintf('%4.2f & %g & %d & %d & %4.2f\\\\% 4.2f\n',0.5*(r'*r),df6,countF6,countN6,tf2,tf1/tf2);
        tfA=[tfA tf/tf2];
        dfA=[dfA df1/df6];
        dim=[dim [num2str(m) 'X' num2str(n)]];
    end
end
semilogy(1:size(dfA,2),dfA);