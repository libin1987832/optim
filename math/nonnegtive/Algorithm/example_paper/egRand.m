% paper rand [1,1] example
%<<<<<<< HEAD
clc
clear
addpath('../FM')
tfA1=[];
dfA1=[];
tfA6=[];
dfA6=[];
dim={};
maIter=500;
for ratio=0.1:0.2:1
    for m=100:100:1000
        
        n=ceil(ratio*m);
        A=2*rand(m,n)-1;
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
        t=clock;
        [xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b,maIter);
        tf1=etime(clock,t);
        %         [xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);
        t=clock;
        [xk6,fk1,xkArr1,countF6,countN6]=hybrid6(x0,A,b,5,10,maIter);
        tf6=etime(clock,t);
        %         [xk7,fk1,xkArr1,countF7,countN7]=hybrid7(x0,A,b,10);
        r=b-A*xk1;
        r(r<0)=0;
        df1=norm(A'*r);
        fprintf('$ %d \\times %d $ & %4.2f & %g & %d & %d & %4.2f &',m,n,0.5*(r'*r),df1,countF1,countN1,tf1);
        r=b-A*xk6;
        r(r<0)=0;
        df6=norm(A'*r);
        fprintf('%4.2f & %g & %d & %d & %4.2f\\\\%% %4.2f %4.2f\n',0.5*(r'*r),df6,countF6,countN6,tf6,tf1/tf6,log10(df1/df6));
        tfA1=[tfA1 tf1];
        tfA6=[tfA6 tf6];
        dfA1=[dfA1 df1];
        dfA6=[dfA6 df6];
        dim=[dim [num2str(m) 'X' num2str(n)]];
    end
end
xais=1:size(tfA1,2);
subplot(1,2,1)
semilogy(1:size(dfA1,2),dfA1,'--o',1:size(dfA1,2),dfA6,'-o');
set(gca,'XTick',xais(1:10:end));
set(gca,'XTicklabel',dim(1:10:end));% X轴的记号
legend("Dax","Ours")
subplot(1,2,2)
plot(1:size(tfA1,2),tfA1,'--o',1:size(tfA6,2),tfA6,'-o');
set(gca,'XTick',xais(1:2:end));
set(gca,'XTicklabel',dim(1:2:end));% X轴的记号
legend("Dax","Ours")
