%% SPL�����˵��㷨��
function [xk,i,h,lamdab]=SPL(A,B,epsilon,epsilon2,e,n,e0)
 x0=zeros(n,1);
 x0(1)=1;
 tic %��ʱ
 i=0;
 options = optimoptions('quadprog','Display','off');%Ϊ�˲���� quadprog�м�Ĺ���
 while 1
lamdab=(x0'*A*x0)/(x0'*B*x0);
yk=-(lamdab*B*x0) ;
xk=BBP(A,yk,e,n,epsilon2);
% xk=quadprog(A,yk,[],[],e,1,e0,[],[],options); %�������ú���quadprog���xk.
dk=xk-x0;
w=(A-lamdab*B)*xk;
gc=(xk'*w)/2;
if abs(gc)<epsilon
break
end     
x0=xk; %��������
i=i+1;
 end 
 toc %��ʱ����
w=(A-lamdab*B)*xk;
h1=min(xk,w);
h=norm(h1); % ������
end