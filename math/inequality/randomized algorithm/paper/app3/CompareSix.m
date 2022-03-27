m1=600;
m2=450;
n=80;
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
xst=zeros(1,4);
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

xs=[x1 x2 x3 x4];
[alpha,beta]
str=['P','R','P','C'];
for i = 1:4
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


