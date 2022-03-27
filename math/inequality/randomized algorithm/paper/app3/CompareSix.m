m1=60;
m2=30;
n=100;
a=0.3;
A1=rand(m1,n);
A2=rand(m2,n);
x=rand(x,0);
d1=A1*x;
d2=A2*x;
dn=floor(m1*a);
u=d1;
u(1:dn,:)=u(1:dn,:)+rand(1:dn,:);

x0=zeros(n,1);

maxIter=20;
iter = 800;
nf = 20;
xst=zeros(1,4);
t=clock;
[x1,arr]=project(L,x2,x0,iter,delt,0,xu,0);
xst(1)=etime(clock,t);
t=clock;
[x2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A,b,x0,maxIter,nf,'RHA');
xst(2)=etime(clock,t);
t=clock;
[x3,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A,b,x0,maxIter,nf,'PHA');
xst(3)=etime(clock,t);
t=clock;
[x4,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A,b,x0,maxIter,nf,'CHA');
xst(4)=etime(clock,t);

xs=[x1 x2 x3 x4];

SNRdB = @(s,n)( 10*log10(sum(s(:).^2)/sum((n(:)-s(:)).^2)) ); 
str=['P','R','P','C'];
for i = 1:4
r = b - A * xs(:,i);
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %sHA &%g &%g & %g & %g \\\\\n', str(i), SNRdB(x,xs(:,i)),r_GS, g_GS, xst(i));
end


