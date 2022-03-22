x=zeros(1,64);
x(1,6:17) = 0.5:0.5:6;
x(1,23) = 9;
x(1,28) = 7;
x(1,38) = 2.1;
x(1,45:55) = ones(1,11)*3.8;
N = length(x);%求取抽样点数
t = (0:N-1);%显示实际时间
hold on
plot(t,x,'g');%绘制时域波形
gausFilter = fspecial('gaussian',[1 3],3);
v1 = ones(1,63)*gausFilter(1);
L1 = diag(v1,1);
v2 = ones(1,64)*gausFilter(2);
L2 = diag(v2);
v3 = ones(1,63)*gausFilter(3);
L3 = diag(v3,-1);
L = L1 + L2 + L3;
u = rand(1,64)*0.15;
x2 = L*x'+ u';
plot(t,x2,'r');
A=[-L;L;eye(64)];
b=[-x2-0.15;x2-0.15;zeros(64,1)];
x0=zeros(64,1);
maxIter = 30;
  [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
  plot(t,xkh','b');
% xff = fft(x);
% stem(1:64, xff);  
  
  % 
% [x,Fs] = wavread('test1');%读取音频数据
% x = x(:,1);
% x = x';
% N = length(x);%求取抽样点数
% t = (0:N-1)/Fs;%显示实际时间
% y = fft(x);%对信号进行傅里叶变换
% f = Fs/N*(0:round(N/2)-1);%显示实际频点的一半
% subplot(211);
% plot(t,x,'g');%绘制时域波形
% axis([0 max(t) -1 1]);
% xlabel('Time / (s)');ylabel('Amplitude');
% title('信号的波形');
% grid;
% subplot(212);
% plot(f,abs(y(1:round(N/2))));
% xlabel('Frequency / (s)');ylabel('Amplitude');
% title('信号的频谱');
% grid;
% %% 以实际频点显示或者显示频率的一半
% f = Fs/N*(0:N-1);
% f = Fs/N*(0:round(N/2)-1);
% %%以实际频点显示，显示的是Fs的一半
% [h,w] = freqz(a,b,512,Fs);
% %% w介于0--pi之间，除以pi，归一化到0--1，乘以奈奎斯特频率，实现显示频率的一半
% [h,w] = freqz(a,b,512);
% f = Fs*w/(2*pi);
% %% 显示实际的时间，抽样点数除以抽样频率
% t = (0:N-1)/Fs;
% %% 信号频率，采样频率，关系
% T/deltat = Fs/fs = T * Fs = N %抽样点数等于信号周期T乘以抽样频率，即每秒抽样多少次，乘以时间