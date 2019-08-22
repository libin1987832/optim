addpath(genpath(pwd));

A=[2,1;-3,-1];
b=[5;-3];
n=2;
e=0.01;

% [x,f]=inexact(b,A,ones(n,1)*100,0.2,0.1,e);

% % % number test
m1=10;m2=600;n=3;density=1;cond=10;delt=1e-5;e=1e-4;
A1=sprand(m1,n,density,1/cond); 
b=sprand(m1,1,density,1/cond)-0.2;
%load matlab.mat
iter=3;
% load matlab_test1.mat
% A1=[1 1;-1 -1];
% b=[1;2]; 
% m1=2;
% y=b-A1*[-sqrt(2)/4;-sqrt(2)/4];
% A1'*y;
% y(y<0)=0;
% gg=0.5*y'*y;

fb=b;
fb(fb<0)=0;
%��ֵ
x0=zeros(n,1)+1;
%��ֵ����� �����Ʋ��Ż��Ľ���ǲ��Ǳ����ֵ��ҪС
init=0.5*fb'*fb
%Q R �ֽ� ���ڹ̶������㷨������ �����ظ��ֽ�
[Q,R]=qr(A1'*A1);
%��¼ ���ַ����ĵĳ�ֵ
x10=x0;
x20=x0;
x30=x0;
%��¼ÿ�ε������Ŀ�꺯��ֵ
fc1=[];
fc2=[];
fc3=[];
%��¼ÿ�ε�������������
nec1=[];
nec2=[];
nec3=[];
for k=1:20
%��¼ÿ�ε�������������
y=b-A1*x10;
All=1:m1;
NE=All(y>=0);
All(NE)=0;
nec1=[nec1;All];
%ǰ����Щ������Ҫ�������������������
[x11,f1]=FM(x10,Q,R,A1,b);
%��¼����Ŀ�꺯��ֵ
fc1=[fc1,f1];

%��¼ÿ�ε�������������
y=b-A1*x20;
All=1:m1;
NE=All(y>=0);
All(NE)=0;
nec2=[nec2;All];
%ǰ����Щ������Ҫ�������������������
[x21,f2]=IFM(x20,A1,b,iter);
%��¼����Ŀ�꺯��ֵ
fc2=[fc2,f2];



y=b-A1*x30;
All=1:m1;
NE=All(y>=0);
All(NE)=0;
nec3=[nec3;All];
[x31,f3]=PAD(x30,A1,b,NE,iter);
% y=b-A1*x1;
% A1'*y
% y(y<0)=0;
% gg2=0.5*y'*y
% NE=All(y>=0)
fc3=[fc3,f3];
% ���µ���ֵ
x10=x11;
x20=x21;
x30=x31;
end
% �ֲ���
x10
x20
x30
% �ֲ���
norm(b-A1*x10)
norm(b-A1*x20)
norm(b-A1*x30)
fcc=[fc1;fc2;fc3]
nec1(end,:)
nec2(end,:)
nec3(end,:)
% A2=sprand(m2,n,density,1/cond);


% b1=rand(m1,1);b2=rand(m2,1);
% A=[A1;-A2];
% b=[b1;-b2];
% % �����������Ӧ�ĺ���ֵ
% fqf(b,A,zeros(n,1))
% delt=0.00001;% delt��ֹ��̬����
% e=0.01; %��ֹ���� �ݶȵķ���С�����ֵ����ֹ
% 
% % ��ȷ������
% %   [x0,f0]=alg1(A,b,ones(n,1)*100,e);
% 
% % �Ǿ�ȷ������
% % inexact(b,A,ones(n,1)*100,0.2,0.1,e); 
% 
%  M=100;% �ͷ�����
% % ��������
%   [x,fk1]=GNP(ones(n,1)*100,M,delt,e,-1*A,-1*b);
