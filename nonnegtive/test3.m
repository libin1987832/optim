addpath(genpath(pwd));

A1=[1,-1;-2,1;1,1;1,1];
b=[2;-1;-2;3];
n=2; 
m1=4;
e=0.01;

% % % number test
% m1=10;m2=600;n=3;density=1;cond=10;delt=1e-5;e=1e-4;
% A1=sprand(m1,n,density,1/cond); 
% b=sprand(m1,1,density,1/cond)-0.2;

iter=3;

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
for k=1:41
%��¼ÿ�ε�������������
y=b-A1*x10;
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec1=[nec1;All];
%ǰ����Щ������Ҫ�������������������
[x11,f1]=FM(x10,Q,R,A1,b);
%��¼����Ŀ�꺯��ֵ
fc1=[fc1,f1];

%��¼ÿ�ε�������������
y=b-A1*x20;
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec2=[nec2;All];
%ǰ����Щ������Ҫ�������������������
[x21,f2]=IFM(x20,A1,b,iter);
%��¼����Ŀ�꺯��ֵ
fc2=[fc2,f2];

y=b-A1*x30;
All=1:m1;
NE=All(y>0.0001);
All(NE)=0;
nec3=[nec3;All];
[x31,f3]=PAD(x30,A1,b,NE,iter);

fc3=[fc3,f3];
% ���µ���ֵ
x10=x11;
x20=x21;
x30=x31;
end
% �ֲ���
disp "����N�κ󣬹̶�������ȷ���������ַ�������ֵ"
xn=[x10,x20,x30]
% �ֲ���
disp "����N�κ󣬹̶�������ȷ���������ַ������պ���ֵ"
% [norm(b-A1*x10),norm(b-A1*x20),norm(b-A1*x30)]
fcc=[fc1;fc2;fc3]
disp "������(0)"
[nec1(end,:);nec2(end,:);nec3(end,:)]
nec3;
fb=b-A1*x30;
fb(fb<0)=0;
A1'*fb
