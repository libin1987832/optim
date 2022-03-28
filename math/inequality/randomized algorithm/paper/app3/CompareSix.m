clear
clc
m1=1500;
m2=1000;
n=250;
al=0.2;
be=0.3;
A1=rand(m1,n);
A2=rand(m2,n);
x=rand(n,1);
d1=A1*x;
d2=A2*x;
dn=floor(m1*al);
u=d1;
u(1:dn,:)=u(1:dn,:)-be*rand(dn,1).*u(1:dn,:);

l=d2;
H=[-A1;A2;eye(n)];
b=[-u;l;zeros(n,1)];
x0=zeros(n,1);

maxIter=100;
iter = 800;
nf = 20;
xst=zeros(1,5);
t=clock;
[x1,alpha,beta,flag]=ablp(A1,A2,u,d2,0.5,0.5,0.1,0.1);
xst(1)=etime(clock,t);
t=clock;
[x2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(H,b,x0,maxIter,nf,'RHA');
xst(2)=etime(clock,t);
t=clock;
[x3,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(H,b,x0,maxIter,nf,'PHA');
xst(3)=etime(clock,t);
t=clock;
[x4,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(H,b,x0,maxIter,nf,'CHA');
xst(4)=etime(clock,t);
maxIter=100;
t=clock;
[x5]=Fluence(A1,A2,u,l,dn,10,maxIter,'quadprog');
xst(5)=etime(clock,t);

xs=[x1 x2 x3 x4 x5];
[alpha,beta]
str=['P','R','P','C','N'];
for i = 1:5
r = -u + A1 * xs(:,i);
rnum1=sum(r>0);
r = l - A2 * xs(:,i);
rnum2=sum(r>0);
r=b-H*xs(:,i);
r(r<0)=0;
r_GS = norm(r);
rum3=sum(x<0);
fprintf('& %sHA &%d &%d &%d & %g & %g \\\\\n', str(i), rnum1,rnum2,rum3,r_GS, xst(i));
end


% ans =
% 
%     0.1000    0.5000
% 
% & PHA &393 &188 &0 & 39.8864 & 7.995 \\
% & RHA &174 &390 &0 & 34.1752 & 0.662 \\
% & PHA &174 &390 &0 & 34.1752 & 0.657 \\
% & CHA &174 &390 &0 & 34.1752 & 0.563 \\
% & NHA &174 &398 &0 & 35.1098 & 0.02 \\
% 
% ans =
% 
%     0.1000    0.4000
% 
% & PHA &802 &56 &0 & 94.9128 & 63.801 \\
% & RHA &316 &664 &0 & 78.9138 & 2.345 \\
% & PHA &316 &664 &0 & 78.9138 & 2.471 \\
% & CHA &316 &664 &0 & 78.9138 & 1.897 \\
% & NHA &292 &688 &0 & 81.3052 & 1.625 \\