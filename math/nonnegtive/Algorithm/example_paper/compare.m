addpath('FM');
addpath('util');
addpath('hybrid');
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
x0=[-1;-0.2];
[Q,R]=qr(A);
xx1=[];
yy1=[];
error=[];
for i=0:1
fk=b-A*x0;
fk(fk<0)=0;
error=[error;norm(A'*(fk))];
xx1=[xx1;x0(1)];
yy1=[yy1;x0(2)];
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
x0=xk;
end
xxF=xx1;
yyF=yy1;
errorF=error;
Qn=Q(:,1:2);
QQ=3*Qn*Qn';
I=diag(ones(4,1))
ssign=getBn(QQ,fm,I);

[xkN,rkN,fk,f0,lambe]=ssqr(x0,A,b);
errorN=[error;norm(A'*rkN)];
% xx1=[xx1;xk(1)];
% yy1=[yy1;xk(2)];
for i=0:3
fk=b-A*x0;
fk(fk<0)=0;
errorF=[errorF;norm(A'*(fk))];
xxF=[xxF;xk(1)];
yyF=[yyF;xk(2)];
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
x0=xk;
end

% plot(4:size(errorF,1),errorF(4:end,1),'-b*',4:size(errorN,1),errorN(4:end),'-or');
d=lineData(A,b,[-2.5,0],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')
% line([-2/3,1],[2/3,-1],'LineStyle','--');
line([-2/3,-2.5],[2/3,2.5],'LineStyle','--');
hold on 
line(xxF',yyF','LineWidth',2)
plot(xxF,yyF,'r*');
H=line([xx1;xkN(1)]',[yy1;xkN(2)]','Color',[.1 .1 .1],'LineWidth',2);
plot([xx1;xkN(1)]',[yy1;xkN(2)]','bo');
line([xx1;]',[yy1;]','Color',[.5 0 0],'LineWidth',2);
    text(xkN(1),xkN(2),'3')
for i=1:size(xxF,1)
    x0=[xxF(i),yyF(i)];
    c=num2str(i);
%     plot(x0(1),x0(2),'*')
    c=[' ',c];
    text(x0(1),x0(2),c)
    %plot(xx1,yy1,'r*');
    hold on
end

[xx1';yy1']
% [xk1,fk1,xkArr1,countF1,countN1]=hybrid6(x0,A,b,3);
% [xk1,fk1,xkArr1,countF1,countN1]=hybrid2(x0,A,b);
% [xk2,fk2,xkArr2,countF2,countN2]=hybrid1(x0,A,b);