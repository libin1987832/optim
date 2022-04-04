 [A1,b1,A2,b2,At1, At2, fm1, fm2, AL1, bL1, AL2, bL2] = readBreast(gamm1);
% [A1,b1,A2,b2,At1, At2, fm1, fm2, AL1, bL1, AL2, bL2] = readHeart(gamm1);
A=A1;
b=b1;
[m,n]=size(A);
x0=zeros(n,1);

maxIter=6000;
iter = 30;
nf = 10;
tol = -1;
xst=zeros(1,4);
t=clock;
xst(1)=etime(clock,t);
t=clock;
[x2,iter,error_k,iter_k,index_k] = GuassSeidelNE(A, b, x0,2 ,maxIter,tol,[],0);
iter
xst(2)=etime(clock,t);
t=clock;
[x3,iter,error_k,iter_k,index_k] = simpleGuassSeidelNE(A, b, x0,2 ,maxIter,tol,[],0);
iter
xst(3)=etime(clock,t);
t=clock;
[x4,iter,error_k,iter_k,index_k] = randGuassSeidelNE(A, b, x0,2 ,maxIter,tol,[],0);
iter
xst(4)=etime(clock,t);

xs=[ x2 x2 x3 x4];

SNRdB = @(s,n)( 10*log10(sum(s(:).^2)/sum((n(:)-s(:)).^2)) ); 
str=['c','C','U','R'];
for i = 2:4
r = b - A * xs(:,i);
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %sCD &%g  & %g & %g \\\\\n', str(i),r_GS, g_GS, xst(i));
end


% & PCD &17.9687 &0.182894 & 0.0986153 & 6.187 \\
% & CCD &17.0119 &0.524659 & 0.082917 & 5.673 \\
% & UCD &4.62 &1.60415 & 0.262694 & 10.364 \\
% & RCD &4.48563 &1.49949 & 0.231798 & 11.358 \\
% 
% & PCD &16.9358 &0.6605 & 0.47467 & 79.913 \\
% & CCD &16.6472 &0.675647 & 0.11328 & 57.519 \\
% & UCD &4.53101 &1.52782 & 0.287233 & 129.173 \\
% & RCD &4.30754 &1.96599 & 0.330704 & 142.153 \\

% 
% & PCD &16.907 &0.667337 & 0.513625 & 78.914 \\
% & CCD &15.5803 &1.28441 & 0.316829 & 29.142 \\
% & UCD &3.98662 &3.13809 & 0.826979 & 65.938 \\
% & RCD &3.96707 &3.78204 & 0.987099 & 72.053 \\