% x0 initial value lambal 
% xk out value 
% r0=(b-Ax0)+ rk=(b-Axk)+ 
% fm=0.5*norm(A*uk-r)^2 fr=norm(yk1)^2 
% fr Ŀǰ�Ĳв� fmͶӰ�ļ����� fk() ���º��Ŀ�꺯�� rk �Ǹ��º������ r0 ���������
function [xk,r0,rk,fk,fm,fr]=BFM(x0,lambal,A,b)
r=b-A*x0;
r(r<0)=0;
r0=r;
fr=0.5*r'*r;
% compute min increase
% B��ֵ
uk=diag(ones(size(x0))*(1/lambal))*(A'*r);
mm=A*uk-r;
fm=0.5*mm'*mm;

xk=x0+uk;
fk=b-A*xk;
fk(fk<0)=0;
rk=fk;
fk=0.5*fk'*fk;

% z0=A*x0-b;
% z0(z0<0)=0;
% u0=A*x0-z0
% 
% zk=A*xk-b;
% zk(zk<0)=0;
% uk=A*xk-zk