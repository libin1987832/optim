clear
clc

n = 500;
w = [1,1];
b = 1;
x = rand( n , 2 );
y = (x*w'-b*ones(n,1) >0);
p1 = x(y,:);
p2 = x(~y,:);
% plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');

%[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
%[xkh2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han([x0;0],A2,b2,maxIter);

%source_train = [p1;p2];
%label_train = [ones(size(p1,1),1);-1*ones(size(p2,1),1)];
gamm = 0.01;
fm1 = size(p1,1);
fm2 = size(p2,1);
A1(1:fm1,:) = p1;
A1(fm1+1:fm1+fm2,:) = -p2;
A2 = [[p1 -1*ones(fm1,1)];[-p2 -ones(fm2,1)]];
b1 = [(1 + gamm)*ones(fm1,1);(-gamm-1)*ones(fm2,1)];% 
b2 = [ones(fm1,1);-ones(fm2,1)];
source_train = [A1(1:fm1,:);-A1(fm1+1:fm1+fm2,:)];
 label_train = [ones(fm1,1);-1*ones(fm2,1)];
%label_train = round(rand(fm1+fm2,1));
SVMModel = fitcsvm(source_train,label_train,'BoxConstraint',30);
wsvn = SVMModel.Beta;
bsvn = SVMModel.Bias;

x0=zeros(size(A1,2),1);
maxIter = 10;
[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
norm(rkh)
[xkh2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han([x0;0],A2,b2,maxIter);
norm(rkh)
ww = [w; wsvn'; xkh1'; ];
bb = [b; -bsvn; 1 ; ];
% 画多条线
d = lineData(ww , bb, [0,1], [0,1]);
plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');
hold on 
line(d(:,[1,2])',d(:,[3,4])')
%line([1,(-bsvn-wsvn(1))/wsvn(2)],[(-bsvn-wsvn(2))/wsvn(1),1])
%plot(d1)
[w/norm(w);wsvn'/norm(wsvn);xkh1'/norm(xkh1);xkh2([1,2])'/norm(xkh2([1,2]))]






% 
%  
% 数据准备
% X = xlsread('C:\Users\user01\Desktop\test.xlsx');
% 二分类 随机生成数据。  200个数据  每个数据2个特征
% data=x;
% label=y;
% data=1*rand(300,2);
% label=zeros(300,1);
% label(sqrt(data(:,1).^2+data(:,2).^2)<8)=1;
% label((data(:,2)+data(:,1)>1))=1;
% 在data上加常数特征项；
% data=[data,ones(size(data,1),1)];
%  
% 打乱循序
% randIndex = randperm(size(data,1));
% data_new=data(randIndex,:);
% label_new=label(randIndex,:);
%  
% 80%训练  20%测试
% k=0.8*size(data,1);
% X1=data_new(1:k,:);
% Y1=label_new(1:k,:);
% X2=data_new(k+1:end,:);
% Y2=label_new(k+1:end,:);
%  
% [m1,n1] = size(X1);
% [m2,n2] = size(X2);
% Features=size(data,2); %特征个数
%  开始训练
% 设定学习率为0.01
% delta=1; 
% lamda=0.2; %正则项系数
%  
% theta1=rand(1,Features);
% theta1=[.5,.5];
% 训练模型
%  
% 梯度下降算法求解theta（每次都是对全部的数据进行训练）
% num = 500; %最大迭代次数
% L=[];
% while(num)
%     dt=zeros(1,Features);
%     loss=0;
%     for i=1:m1
%         xx=X1(i,1:Features);
%         yy=Y1(i,1);
%         h=1/(1+exp(-(theta1 * xx')));
%         dt=dt+(h-yy) * xx;
%         loss=loss+ yy*log(h)+(1-yy)*log(1-h);
%     end
%     loss=-loss/m1;
%     L=[L,loss];
%      
%     theta2=theta1 - delta*dt/m1 - lamda*theta1/m1;
%     theta1=theta2;
%     num = num - 1;
%      
%     if loss<0.01
%         break;
%     end
% end
% figure
% subplot(1,2,1)
% plot(L)
% title('loss')
% theta1 
% subplot(1,2,2)
% x=0:0.01:1;
% y=(-theta1(1)*x-theta1(3))/theta1(2);
% plot(x,y,'linewidth',2)
% hold on
% plot(data(label==1,1),data(label==1,2),'ro')
% hold on
% plot(data(label==0,1),data(label==0,2),'go')
% axis([0 1 0 1])
%  
 
%测试数据
% acc=0;
% for i=1:m2
%     xx=X2(i,1:Features)';
%     yy=Y2(i);
%     finil=1/(1+exp(-theta2 * xx));
%     if finil>0.5 && yy==1
%         acc=acc+1;
%     end
%     if finil<=0.5 && yy==0
%         acc=acc+1;
%     end
% end
% acc/m2