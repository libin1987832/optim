clear
clc

n = 200;
w = rand(1,2);
b = rand(1);
x = rand( n , 2 );
y = (x*w'-b*ones(n,1) >0);
p1 = x(y,:);
p2 = x(~y,:);
% plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');

%[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
%[xkh2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han([x0;0],A2,b2,maxIter);

source_train = [p1;p2];
label_train = [ones(size(p1,1),1);-1*ones(size(p2,1),1)];
%label_train = round(rand(fm1+fm2,1));
SVMModel = fitcsvm(source_train,label_train,'BoxConstraint',30);
wsvn = SVMModel.Beta;
bsvn = SVMModel.Bias;
% ª≠∂‡Ãıœﬂ
d = lineData(wsvn' , -bsvn, [0,1], [0,1]);
plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');
hold on 
line(d(1,[1,2])',d(1,[3,4])')
%line([1,(-bsvn-wsvn(2))/wsvn(1)],[(-bsvn-wsvn(1))/wsvn(2),1])
%plot(d1)
