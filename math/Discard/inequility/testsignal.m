x=zeros(1,64);
x(1,6:17) = 0.5:0.5:6;
x(1,23) = 9;
x(1,28) = 7;
x(1,38) = 2.1;
x(1,45:55) = ones(1,11)*3.8;
N = length(x);%��ȡ��������
t = (0:N-1);%��ʾʵ��ʱ��
hold on
plot(t,x,'g');%����ʱ����
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
% [x,Fs] = wavread('test1');%��ȡ��Ƶ����
% x = x(:,1);
% x = x';
% N = length(x);%��ȡ��������
% t = (0:N-1)/Fs;%��ʾʵ��ʱ��
% y = fft(x);%���źŽ��и���Ҷ�任
% f = Fs/N*(0:round(N/2)-1);%��ʾʵ��Ƶ���һ��
% subplot(211);
% plot(t,x,'g');%����ʱ����
% axis([0 max(t) -1 1]);
% xlabel('Time / (s)');ylabel('Amplitude');
% title('�źŵĲ���');
% grid;
% subplot(212);
% plot(f,abs(y(1:round(N/2))));
% xlabel('Frequency / (s)');ylabel('Amplitude');
% title('�źŵ�Ƶ��');
% grid;
% %% ��ʵ��Ƶ����ʾ������ʾƵ�ʵ�һ��
% f = Fs/N*(0:N-1);
% f = Fs/N*(0:round(N/2)-1);
% %%��ʵ��Ƶ����ʾ����ʾ����Fs��һ��
% [h,w] = freqz(a,b,512,Fs);
% %% w����0--pi֮�䣬����pi����һ����0--1�������ο�˹��Ƶ�ʣ�ʵ����ʾƵ�ʵ�һ��
% [h,w] = freqz(a,b,512);
% f = Fs*w/(2*pi);
% %% ��ʾʵ�ʵ�ʱ�䣬�����������Գ���Ƶ��
% t = (0:N-1)/Fs;
% %% �ź�Ƶ�ʣ�����Ƶ�ʣ���ϵ
% T/deltat = Fs/fs = T * Fs = N %�������������ź�����T���Գ���Ƶ�ʣ���ÿ��������ٴΣ�����ʱ��