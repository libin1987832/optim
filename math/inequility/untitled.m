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
output = zeros(2,6);
for i = 1:2
    count1 = sum(At1*xkh(:,i) < gamm(1,i));
    count2 = size(At1,1);
%     percent1 = count1/count2;
%     [-count1+count2 percent1]
    count21 = sum(At2*xkh(:,i) > gamm(1,i));
    count22 = size(At2,1);
%     percent2 = count21/count22;
%     [-count21+count22 percent2]
    sumerror = -count21+count22-count1+count2;
    sumcount = count2+count22;
    output(1,3*i-2:i*3) = [sumerror, sumcount, 1-sumerror/sumcount];
end

AL1 = [-A1,-eye(size(A1,1))];
bL1 = -b1;
fL1 = [zeros(1, 9), ones(1, fm1)/fm1, ones(1, fm2)/fm2];
lbL1 = [-ones(9,1)*Inf; zeros(fm1+fm2 , 1)];
[xkh1, f1, exit1] = linprog(fL1,AL1,bL1,[],[],lbL1,[],x0);

AL2 = [-A2,-eye(size(A2,1))];
bL2 = -b2;
fL2 = [zeros(1, 10), ones(1, fm1)/fm1, ones(1, fm2)/fm2];
lbL2 = [-ones(10,1)*Inf; zeros(fm1+fm2 , 1)];
[xkh2, f2, exit2] = linprog(fL2,AL2,bL2,[],[],lbL2,[],[x0;0]);
xw = [xkh; gamm]; 
xkh = [xkh1(1:9,1) xkh2(1:9,1)];
gamm = [gamm1,xkh2(10,1)];
xw = [xw [xkh;gamm]];
for i = 1:2
    count1 = sum(At1*xkh(:,i) < gamm(1,i));
    count2 = size(At1,1);
    percent1 = count1/count2;
%     [-count1+count2 percent1]
    count21 = sum(At2*xkh(:,i) > gamm(1,i));
    count22 = size(At2,1);
%     percent2 = count21/count22;
    sumerror = -count21+count22-count1+count2;
    sumcount = count2+count22;
    output(2,3*i-2:i*3) = [sumerror, sumcount, 1-sumerror/sumcount];
end


source_train = [-1*A1(1:fm1,:);A1(fm1+1:fm1+fm2,:)];
 label_train = [ones(fm1,1);-1*ones(fm2,1)];
%label_train = round(rand(fm1+fm2,1));
SVMModel = fitcsvm(source_train,label_train,'BoxConstraint',30);
[ans1,~]=predict(SVMModel, At1);
[ans2,~]=predict(SVMModel, At2);
errorcount1 = sum(ans1 < 0.5);
errorcount2 = sum(ans2 > 0.5);
sumerror = errorcount1 + errorcount2;
output
output2 = [sumerror,sumcount,1-sumerror/sumcount]
xw = [xw [SVMModel.Beta;SVMModel.Bias]]'
beta = xw(:,1:end-1);
norms=repmat(sqrt(sum(beta.^2,2)),1,9);
beta = beta./norms;
beta*beta'

% m = size(A1,1);
% k = size()
% f=[zeros(1,size(xkh1, 1)) 0 ones(size(A1,))]




%flops(0)           %start global flop count at 0
% A=[1 2 3; 4 5 6];
% b=[7 8 9]';
% x=A*b;             %do the operation
% addflops(A*b) %do the counting
%flops              %display count so far
