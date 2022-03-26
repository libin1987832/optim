%%%%1
x=zeros(1,64);
x(1,6:17) = 0.5:0.5:6;
x(1,23) = 9;
x(1,28) = 7;
x(1,38) = 2.1;
x(1,45:55) = ones(1,11)*3.8;
N = length(x);%求取抽样点数
t = (0:N-1);%显示实际时间

gausFilter = fspecial('gaussian',[1 3],3);
v1 = ones(1,63)*gausFilter(1);
L1 = diag(v1,1);
v2 = ones(1,64)*gausFilter(2);
L2 = diag(v2);
v3 = ones(1,63)*gausFilter(3);
L3 = diag(v3,-1);
L = L1 + L2 + L3;
delt=0.15;
u = rand(1,64)*delt;
x2 = L*x'+ u';
xu=12;
A=[-L;L;eye(64);-eye(64)];
b=[-x2-0.15;x2-0.15;zeros(64,1);-ones(64,1)*12];
x0=zeros(64,1);

maxIter=100;
iter = 1000;
nf = 5;
xst=zeros(1,4);
t=clock;
[x1,arr]=project(L,b,x0,iter,delt,0,xu,0);
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
str=["Proj","RHA",'PHA','CHA'];
for i = 1:4
r = b - A * xs(:,i);
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %s & %g & %g & %g \\\\\n', cell2str(str(i)), r_GS, g_GS, xst(i));
end


% % SNRdB = @(s,n)( 10*log10((sum((n(:).^2)-sum((s(:).^2))))/sum(s(:).^2))); 
% SNRdB(x,x2)
% SNRdB(x,xkh)
% SNRdB(x,xpp)


