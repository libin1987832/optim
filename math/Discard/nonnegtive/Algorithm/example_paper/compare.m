% compare our hybrid and Dax
% draw picture to compare the different methods
addpath('../FM');
addpath('../util');
addpath('../hybrid');
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
x0=[-1;-0.2];
[Q,R]=qr(A);
xx1=[x0(1)];
yy1=[x0(2)];
ssign=[];
skn=[];
fk=b-A*x0;

Qn=Q(:,1:2);
QQ=2*Qn*Qn';
I=diag(ones(4,1))
ssign1=getBn(QQ,fk,I);
ssign=[ssign1];
rkn=getEstimatedR(QQ,fk,I);
skn=[rkn];
fk(fk<0)=0;
fkk1=[0.5*(fk'*fk)];
error=[norm(A'*(fk))];

% matlab have four column

% take two steps for both algorithm
for i=1:1
fk=b-A*x0;
fk(fk<0)=0;
fkk1=[fkk1;(fk'*fk)];
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
x0=xk;
xx1=[xx1;x0(1)];
yy1=[yy1;x0(2)];
error=[error;norm(A'*(rk))];
end
xxF=xx1;
yyF=yy1;
errorF=error;

% predict
ssign1=getBn(QQ,fm,I);
ssign=[ssign ssign1];
rkn=getEstimatedR(QQ,fm,I);
skn=[skn rkn];

% take newton type method
[xkN,fk,dh,rkk]=ssqr4(x0,A,b);
errorN=[error;norm(A'*rkk)];
xx1=[xx1;xkN(1)];
yy1=[yy1;xkN(2)];
for i=0:3
fk=b-A*x0;
fk(fk<0)=0;
fkk1=[fkk1;(fk'*fk)];
[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
errorF=[errorF;norm(A'*(rk))];
xxF=[xxF;xk(1)];
yyF=[yyF;xk(2)];
x0=xk;
end

% plot(4:size(errorF,1),errorF(4:end,1),'-b*',4:size(errorN,1),errorN(4:end),'-or');
% plot problem
d=lineData(A,b,[-2.5,0],[-1.5,3]);
line(d(:,1:2)',d(:,3:4)')

% line([-2/3,1],[2/3,-1],'LineStyle','--');
% plot solution 
line([-2/3,-2.5],[2/3,2.5],'LineStyle','--');
hold on 
% plot Dax
line(xxF',yyF','LineWidth',2,'linestyle','--','color','b')
plot(xxF,yyF,'b*');
H=line([xx1;xkN(1)]',[yy1;xkN(2)]','Color',[.1 .1 .1],'LineWidth',2,'linestyle','--');

plot([xx1;xkN(1)]',[yy1;xkN(2)]','ko');
line([xx1;]',[yy1;]','Color','k','LineWidth',2);
    text(xkN(1)+0.05,xkN(2),'2','color','k')
line([xx1(1:2)]',[yy1(1:2)]','Color','r','LineWidth',2);
for i=1:size(xxF,1)-2
    x0=[xxF(i),yyF(i)];
    c=num2str(i-1);
%     plot(x0(1),x0(2),'*')
    c=[' ',c];
    if i<3
    text(x0(1)+0.05,x0(2)-0.01,c,'color','r')
    else
    text(x0(1)+0.05,x0(2)-0.01,c,'color','b')   
    end    %plot(xx1,yy1,'r*');
    hold on
end

[xx1';yy1']
skn
for i=0:2
fprintf('$ x_%i $ & (%4.2f,%4.2f) & %4.2f & {1,2}& (%4.2f,%4.2f) & %4.2f & {1,2} & - \\\\\n',i,xxF(i+1),yyF(i+1),errorF(i+1)...
    ,xx1(i+1),yy1(i+1),errorN(i+1));
end
for i=3:5
fprintf('$ x_%i $ & (%4.2f,%4.2f) & %4.2f & {1,2} & - & - & - & -\\\\\\n',i,xxF(i+1),yyF(i+1),errorF(i+1));
end

% [xk1,fk1,xkArr1,countF1,countN1]=hybrid6(x0,A,b,3);
% [xk1,fk1,xkArr1,countF1,countN1]=hybrid2(x0,A,b);
% [xk2,fk2,xkArr2,countF2,countN2]=hybrid1(x0,A,b);