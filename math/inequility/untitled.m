addpath('./util')
gamm1 = 0.5;
[A1,b1,A2,b2,At1, At2, fm1, fm2, AL1, bL1, AL2, bL2] = readBreast(gamm1);
maxIter = 900;
str = ['D','C','R','P'];
x0=zeros(size(A1,2),1);
[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
[xkh2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han([x0;0],A2,b2,maxIter);
xkh = [xkh1 xkh2(1:end-1)];
gamm = [gamm1,xkh2(end,1)];
for i = 1:2
    count1 = sum(At1*xkh(:,i) < gamm(1,i));
    count2 = size(At1,1);
    percent1 = count1/count2;
    [-count1+count2 percent1]
    count21 = sum(At2*xkh(:,i) > gamm(1,i));
    count22 = size(At2,1);
    percent2 = count21/count22;
    [-count21+count22 percent2]
end

AL1 = [-A1,-eye(size(A1,1))];
bL1 = -b1;
fL1 = [zeros(1, 9), ones(1, fm1)/fm1, ones(1, fm2)/fm2];
lbL1 = [-ones(9,1)*Inf, zeros(fm1+fm2 , 1)];

% m = size(A1,1);
% k = size()
% f=[zeros(1,size(xkh1, 1)) 0 ones(size(A1,))]




%flops(0)           %start global flop count at 0
% A=[1 2 3; 4 5 6];
% b=[7 8 9]';
% x=A*b;             %do the operation
% addflops(A*b) %do the counting
%flops              %display count so far
