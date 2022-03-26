%%%%1
% x=zeros(1,640);
% x(1,60:170) = 0.5:0.5:6;
% x(1,230) = 9;
% x(1,280) = 7;
% x(1,380) = 2.1;
% x(1,450:550) = ones(1,11)*3.8;
N=1000;
x=zeros(1,N);
x(1,60:170) = 5:(6-5)/110:6;
x(1,230:300) = 9;
x(1,380:420) = 7;
x(1,430) = 2.1;
x(1,450:550) = ones(1,101)*3.8;
N = length(x);%求取抽样点数
t = (0:N-1);%显示实际时间

% gausFilter = fspecial('gaussian',[1 3],3);
m=20;
 gausFilter = fspecial('gaussian',[1 m],3);
 L=zeros(N,N);
 for i = 1:m
   v1=  ones(1,N-i+1) * gausFilter(i);
   L1 = diag(v1,i-1);
   L= L + L1;
% v1 = ones(1,N-1)*gausFilter(1);
% L1 = diag(v1,1);
% v2 = ones(1,N)*gausFilter(2);
% L2 = diag(v2);
% v3 = ones(1,N-1)*gausFilter(3);
% L3 = diag(v3,-1);
 end
% L = L1 + L2 + L3;
L=L+L';
 delt=0.15;
u=zeros(1,N);
uN=floor(N);
u(1:uN) = rand(1,uN)*0.9*delt;
x2 = L*x'+ u';
xu=12;
A=[-L;L;eye(N);-eye(N)];
b=[-x2-delt;x2-delt;zeros(N,1);-ones(N,1)*xu];
x0=zeros(N,1);

maxIter=30;
iter = 300;
nf = 10;
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


% & PHA &16.6126 &1.5116e-15 & 7.3582e-16 & 0.839 \\
% & RHA &16.908 &4.45095e-16 & 2.56656e-16 & 0.013 \\
% & PHA &16.9011 &6.28876e-16 & 3.63282e-16 & 0.015 \\
% & CHA &16.0298 &6.2865e-16 & 3.63078e-16 & 0.011 \\

% & PHA &20.0875 &9.83359e-15 & 4.2915e-15 & 1.983 \\
% & RHA &19.948 &2.67462e-15 & 9.12252e-16 & 1.237 \\
% & PHA &19.7099 &3.00804e-15 & 1.68034e-15 & 1.283 \\
% & CHA &18.6069 &2.60379e-15 & 1.23702e-15 & 1.218 \\