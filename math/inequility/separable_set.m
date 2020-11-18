clear
clc
gamm = 0.6;
n = 300;
w = [1,1]+rand(1,2);
b = 1+rand(1);
x = rand( n , 2 );
%% triangle
% w1 = [1,1];w2 = [1,-1];
% b1 = 0.2;b2 = 1;
% y1 = (x*w1'-(b1+1)*ones(n,1) < 0 & x*w2'-b1*ones(n,1) < 0);
% y2 = (x*w1'-(1-b1)*ones(n,1) > 0 & x*w2'+b1*ones(n,1) > 0);
% p3 = x(y1&y2,:);
% row = size(p3,1)/2;
% p1 = x(y1&~y2,:);
% p2 = x(~y1&y2,:);
% p1 = [p1;p3(1:row,:)];
% p2 = [p2;p3(row+1:end,:)];

%% guass nosise
% y = (x*w'-b*ones(n,1) >0);
% 
% p1 = x(y,:);
% p2 = x(~y,:);
% p1= p1 + mvnrnd([0,0],[0.005,0;0,0.005],size(p1,1)); 
% p2= p2 + mvnrnd([0,0],[0.005,0;0,0.005],size(p2,1));

%% inconsistent data
%y = (x(:,1)>0.5) | (x(:,2)>0.5); 
%% consistent data 
y = (x*w'-b*ones(n,1) >0);

p1 = x(y,:);
p2 = x(~y,:);

%% test data
% px1 = 0:0.2:1;
% p1 = [px1' -px1'+1.5];
% p2 = [px1' -px1'+0.5];

%%
fm1 = size(p1,1);
fm2 = size(p2,1);
A1(1:fm1,:) = p1;
A1(fm1+1:fm1+fm2,:) = -p2;
A2 = [[p1 -1*ones(fm1,1)];[-p2 ones(fm2,1)]];
% b1 = [(1 + gamm)*ones(fm1,1);(-1-gamm)*ones(fm2,1)];% 
b2 = [ones(fm1,1);ones(fm2,1)];
source_train = [A1(1:fm1,:);-A1(fm1+1:fm1+fm2,:)];
label_train = [ones(fm1,1);-1*ones(fm2,1)];
%label_train = round(rand(fm1+fm2,1));
tic
SVMModel = fitcsvm(source_train,label_train,'BoxConstraint',30);
toc
wsvn = SVMModel.Beta;
bsvn = SVMModel.Bias;

x0=zeros(size(A1,2),1);
maxIter = 10;
% there are some errors;
%[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
tic
[xkh2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han([x0;0],A2,b2,maxIter);
toc

numberOfbeta=2;
AL2 = [-A2,-eye(size(A2,1))];
bL2 = -b2;
fL2 = [zeros(1, numberOfbeta+1), ones(1, fm1)/fm1, ones(1, fm2)/fm2];
lbL2 = [-ones(numberOfbeta+1,1)*Inf; zeros(fm1+fm2 , 1)];
tic
[xkh2, f2, exit2] = linprog(fL2,AL2,bL2,[],[],lbL2,[],[x0;0]);
toc


ww = [w; wsvn'; xkh2([1,2])';xkh2([1,2])'; xkh2([1,2])' ];
bb = [b; -bsvn; xkh2(3)+1 ;xkh2(3)-1;xkh2(3)];
display = [5];
% ��������
d = lineData(ww(display,:) , bb(display,:) , [0,1], [0,1]);
plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');
hold on 
line(d(:,[1,2])',d(:,[3,4])')
%line([1,(-bsvn-wsvn(1))/wsvn(2)],[(-bsvn-wsvn(2))/wsvn(1),1])
[[w b]/norm(w);[wsvn' -bsvn]/norm(wsvn);[xkh2([1,2])' xkh2(3)]/norm(xkh2([1,2]))];
 maxIter = 20;
 nf = 3;
 str = ['D','C','R','P'];
 steplength = 0;%1/(max(eig(ATA))+0.0001);
    for i=1:4
        type = str(i);
        tic
        [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A2,b2,[x0;0],steplength,maxIter,nf,[type,'HA']);
        toc
        %fprintf('& %s & %g \\\\\n',[type,'HA'],tfD);
    end
     
  