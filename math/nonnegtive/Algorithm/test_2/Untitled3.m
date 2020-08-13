clear
f1=load('ff1.mat');
f2=load('ff2.mat');
f3=load('ff3.mat');
f4=load('ff4.mat');
f5=load('ff5.mat');
f6=load('ff6.mat');
f7=load('ff7.mat');
f8=load('ff8.mat');
f9=load('ff9.mat');
f10=load('ff10.mat');
summ=[f1.Arecord(:,1:3),f2.Arecord(:,3),f3.Arecord(:,3),f4.Arecord(:,3),...
    f5.Arecord(:,3),f6.Arecord(:,3),f7.Arecord(:,3),f8.Arecord(:,3)...
    ,f9.Arecord(:,3),f10.Arecord(:,3)];
summ2=[summ mean(summ(:,3:12),2) std(summ(:,3:12),0,2)];
subplot(1,3,1)
plot(0.1:0.1:1,summ2(1:10,13))
subplot(1,3,2)
plot(0.1:0.1:1,summ2(11:20,13))
subplot(1,3,3)
plot(0.1:0.1:1,summ2(21:30,13))