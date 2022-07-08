clc
clear
str = 'bcsstk02.mtx bcsstk04.mtx';
strindex = [1 12; 14 25];
iterstr = [98 422;984 973];
for i = 1:size(strindex,1) 
[B,rows,cols,entries,rep,field,symm] = mmread(str(1,strindex(i,1):strindex(i,2)));
n=cols;
d = sprand(n,1,0.5);
x1=sparse(ones(n,1)./n);
t=clock;
a = (d' * x1) / (x1' * B *x1);
[x, iter, Abpp] =  FLQP(a, B, d, 0, n, 2, 100, 20, 1e-6, 1e-9);
tf=etime(clock,t);
xBx = x' * B * x;
dx = d'* x;
ay = dx / xBx;
fa = dx - ay * xBx;
disp(['flqp:it=' num2str(iter) ',a=' num2str(ay) ',fa='...
    num2str(fa) ',abpp=' num2str(Abpp) ',nLsyst=' num2str(Abpp * iter) ',Cpu=' num2str(tf)])
end