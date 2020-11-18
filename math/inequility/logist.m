
clear
clc
 
%% 数据准备
%X = xlsread('C:\Users\user01\Desktop\test.xlsx');
%二分类 随机生成数据。  200个数据  每个数据2个特征
%data=1*rand(300,2);
%label=zeros(300,1);
%label(sqrt(data(:,1).^2+data(:,2).^2)<8)=1;
%label((data(:,2)+data(:,1)>1))=1;
n = 500;
w = [1,1]+rand(1,2);
b = 1+rand(1);
x = rand( n , 2 );

w1 = [1,1];w2 = [1,-1];
b1 = 0.2;b2 = 1;
y1 = (x*w1'-(b1+1)*ones(n,1) < 0 & x*w2'-b1*ones(n,1) < 0);
y2 = (x*w1'-(1-b1)*ones(n,1) > 0 & x*w2'+b1*ones(n,1) > 0);
p3 = x(y1&y2,:);
row = size(p3,1)/2;
p1 = x(y1&~y2,:);
p2 = x(~y1&y2,:);
p1 = [p1;p3(1:row,:)];
p2 = [p2;p3(row+1:end,:)];
data = [p1;p2];
data=[data,ones(size(data,1),1)];%在data上加常数特征项；
label =[ones(size(p1,1),1);zeros(size(p2,1),1)];
 
randIndex = randperm(size(data,1));%打乱循序

data_new=data(randIndex,:);
label_new=label(randIndex,:);
 
k=0.8*size(data,1);% 80%训练  20%测试

X1=data_new(1:k,:);
Y1=label_new(1:k,:);
X2=data_new(k+1:end,:);
Y2=label_new(k+1:end,:);
 
[m1,n1] = size(X1);
[m2,n2] = size(X2);
Features=size(data,2); %特征个数

delta=1;%设定学习率为0.01 
lamda=0.2; %正则项系数
 
theta1=rand(1,Features);
theta1=[.5,.5,0];
%训练模型
 
num = 300; %最大迭代次数 %梯度下降算法求解theta（每次都是对全部的数据进行训练）

L=[];
while(num)
    dt=zeros(1,Features);
    loss=0;
    for i=1:m1
        xx=X1(i,1:Features);
        yy=Y1(i,1);
        h=1/(1+exp(-(theta1 * xx')));
        dt=dt+(h-yy) * xx;
        loss=loss+ yy*log(h)+(1-yy)*log(1-h);
    end
    loss=-loss/m1;
    L=[L,loss];
     
    theta2=theta1 - delta*dt/m1 - lamda*theta1/m1;
    theta1=theta2;
    num = num - 1;
     
    if loss<0.01
        break;
    end
end
figure
subplot(1,2,1)
plot(L)
title('loss')
 
subplot(1,2,2)
x=0:0.1:10;
y=(-theta1(1)*x-theta1(3))/theta1(2);
plot(x,y,'linewidth',2)
hold on
plot(data(label==1,1),data(label==1,2),'ro')
hold on
plot(data(label==0,1),data(label==0,2),'go')
axis([0 1 0 1])
 
 
acc=0;%测试数据

for i=1:m2
    xx=X2(i,1:Features)';
    yy=Y2(i);
    finil=1/(1+exp(-theta2 * xx));
    if finil>0.5 && yy==1
        acc=acc+1;
    end
    if finil<=0.5 && yy==0
        acc=acc+1;
    end
end
acc/m2
 