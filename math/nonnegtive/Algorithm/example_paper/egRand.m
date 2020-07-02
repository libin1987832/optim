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
 for m=1000:100:1000
for ratio=0.1:0.01:0.3
        n=ceil(ratio*m);
        A=2*rand(m,n)-1;
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
        t=clock;
        [xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b,maIter);
        tf1=etime(clock,t);
        %         [xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);
        t=clock;
        [xk6,fk6,xkArr6,countF6,countN6]=hybrid6(x0,A,b,5,20,maIter);
        [xk61,fk61,xkArr61,countF61,countN61]=hybrid6(x0,A,b,5,20,maIter,1);
%         if n > 60
%         [xk6,fk1,xkArr1,countF6,countN6]=hybrid6(x0,A,b,5,20,maIter);
%         else
%         [xk6,fk1,xkArr1,countF6,countN6]=hybrid6(x0,A,b,30,50,maIter); 
%         end
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
        dim=[dim num2str(n)];
    end
end
xais=1:size(tfA1,2);
itxais=5;
subplot(1,2,1)
%plot(1:size(dfA1,2),dfA1,'--o',1:size(dfA1,2),dfA6,'-o');
semilogy(1:size(dfA1,2),dfA1,'--o',1:size(dfA1,2),dfA6,'-o');
set(gca,'XTick',xais(1:itxais:end));
set(gca,'XTicklabel',dim(1:itxais:end));% X轴的记号
legend("Dax","Ours")
xlabel('varibles');
ylabel('the log of the norm of gradient ');
title('the accuracy as the varibles increase');
subplot(1,2,2)
plot(1:size(tfA1,2),tfA1,'--o',1:size(tfA6,2),tfA6,'-o');
set(gca,'XTick',xais(1:itxais:end));
set(gca,'XTicklabel',dim(1:itxais:end));% X轴的记号
legend("Dax","Ours")
xlabel('varibles');
ylabel('the time(s)');
title('the cost time as the varibles increase');
save('data6','dfA1','dfA6','tfA1','tfA6')