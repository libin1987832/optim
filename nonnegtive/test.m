addpath(genpath(pwd));

% % number test
m1=500;m2=500;n=1000;density=0.1;cond=100;delt=1e-5;e=1e-4;
A1=sprand(m1,n,density,1/cond);
A2=sprand(m2,n,density,1/cond);
b1=rand(m1,1);b2=rand(m2,1);
A=[A1;-A2];
b=[b1;-b2];
% �����������Ӧ�ĺ���ֵ
fqf(b,A,zeros(n,1))
delt=0.00001;% delt��ֹ��̬����
e=0.01; %��ֹ���� �ݶȵķ���С�����ֵ����ֹ

% ��ȷ������
 [x0,f0]=alg1(A,b,ones(n,1)*100,e);

% �Ǿ�ȷ������
%inexact(b,A,ones(n,1)*100,0.2,0.1); 

 M=1000;% �ͷ�����
% ��������
%  [x,fk1]=GNP(ones(n,1)*100,M,delt,e,-1*A,-1*b);
