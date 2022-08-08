%%%%1
x=zeros(1,64);
x(1,6:17) = 0.5:0.5:6;
x(1,23) = 9;
x(1,28) = 7;
x(1,38) = 2.1;
%x(1,45:50) = ones(1,6)*3.8;
x(1,45) = 7.8;
x(1,50) = 8.8;
x(1,55) = 5.8;
N = length(x);%��ȡ��������
t = (0:N-1);%��ʾʵ��ʱ��
hold on
subplot(411)
plot(t,x,'g');%����ʱ����
title('ԭʼ�ź�')
gausFilter = fspecial('gaussian',[1 3],3);
v1 = ones(1,63)*gausFilter(1);
L1 = diag(v1,1);
v2 = ones(1,64)*gausFilter(2);
L2 = diag(v2);
v3 = ones(1,63)*gausFilter(3);
L3 = diag(v3,-1);
L = L1 + L2 + L3;
ur=0.3;
u = rand(1,64)*ur;
x2 = L*x'+ u';
subplot(412)
plot(t,x2,'r');
title('�����ź�')
A=[-L;L;eye(64);-eye(64)];
b=[-x2-ur;x2-ur;zeros(64,1);-ones(64,1)*12];
x0=zeros(64,1);
maxIter = 100;
nf=5;
[xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A,b,x0,maxIter,nf,'PHA');
rr=b-A*xkh;
rr(rr<0)=0;
Ar=A'*rr;
norm(Ar)
rr=b-A*x';
rr(rr<0)=0;
Ar=A'*rr;
norm(Ar)
rr=b-A*x2;
rr(rr<0)=0;
Ar=A'*rr;
norm(Ar)
[xpp,arr]=project(L,x2,x0,1000,ur,0,12,0);
rr=b-A*xpp;
rr(rr<0)=0;
Ar=A'*rr;
norm(Ar)
subplot(413)
plot(t,xkh','b');
title('�ع��ź�')
subplot(414)
plot(t,x,'g',t,x2,'r',t,xkh,'b')
%  plot(t,xpp,'b')
legend('ԭʼ�ź�','�����ź�','�ع��ź�')
title('�����źŶԱ�')
SNRdB = @(s,n)( 10*log10(sum(s(:).^2)/sum((n(:)-s(:)).^2)) ); 
% SNRdB = @(s,n)( 10*log10((sum((n(:).^2)-sum((s(:).^2))))/sum(s(:).^2))); 
% SNRdB(x,x2)
% SNRdB(x,xkh)
% SNRdB(x,xpp)
% SNR_singlech(x,x)


% function y = add_gaussian_noise_snr_db(signal, snr)
% x = signal(:)';
% power_of_signal = (x*x')/length(x);
% power_of_noise = power_of_signal/(10^(snr/10)); % the std of the noise is the power of the noise
% noise0 = randn(size(x));
% noise1 = (noise0 - mean(noise0))/std(noise0 - mean(noise0));
% y = sqrt(power_of_noise)*noise1+ x;
% end

% r=b-A*xkh;
% r(r<0)=0;
% r'*r
% r=b-A*xpp;
% r(r<0)=0;
% r'*r
%  
 %%%%%2
%  x=zeros(1,64);
%  x=0:2*pi/64:2*pi;   %��ʼ����������ֹ��ֵ
% x=x(1:64);
% x=(1+1/2*cos(2*pi*1000*x)).*cos(2*pi*10000*x)*2;%���ƺ���
% t=-2*pi/100:pi/1024:2*pi/100;
% 
% y=square(2*pi*30*t,50)+abs(min(y));
% x=[y y(1:23)];
% N = length(x);%��ȡ��������
% t = (0:N-1);%��ʾʵ��ʱ��
% hold on
% subplot(434)
% plot(t,x,'g');%����ʱ����
% gausFilter = fspecial('gaussian',[1 3],3);
% v1 = ones(1,63)*gausFilter(1);
% L1 = diag(v1,1);
% v2 = ones(1,64)*gausFilter(2);
% L2 = diag(v2);
% v3 = ones(1,63)*gausFilter(3);
% L3 = diag(v3,-1);
% L = L1 + L2 + L3;
% u = rand(1,64)*0.15;
% x2 = L*x'+ u';
% subplot(435)
% plot(t,x2,'r');
% A=[-L;L;eye(64);-eye(64)];
% b=[-x2-0.15;x2-0.15;zeros(64,1);-ones(64,1)*12];
% x0=zeros(64,1);
% maxIter = 30;
% [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
% subplot(436)
%  plot(t,xkh','b');
%  %%%%%3
%   x=zeros(1,64);
% x(1,6:17) = 0.5:0.5:6;
% x(1,23) = 9;
% x(1,28) = 7;
% x(1,38) = 2.1;
% x(1,45:55) = ones(1,11)*3.8;
% N = length(x);%��ȡ��������
% t = (0:N-1);%��ʾʵ��ʱ��
% hold on
% plot(t,x,'g');%����ʱ����
% gausFilter = fspecial('gaussian',[1 3],3);
% v1 = ones(1,63)*gausFilter(1);
% L1 = diag(v1,1);
% v2 = ones(1,64)*gausFilter(2);
% L2 = diag(v2);
% v3 = ones(1,63)*gausFilter(3);
% L3 = diag(v3,-1);
% L = L1 + L2 + L3;
% u = rand(1,64)*0.15;
% x2 = L*x'+ u';
% plot(t,x2,'r');
% A=[-L;L;eye(64);-eye(64)];
% b=[-x2-0.15;x2-0.15;zeros(64,1);-ones(64,1)*12];
% x0=zeros(64,1);
% maxIter = 30;
% [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
%  plot(t,xkh','b');
%  %%%%%%%4
%   x=zeros(1,64);
% x(1,6:17) = 0.5:0.5:6;
% x(1,23) = 9;
% x(1,28) = 7;
% x(1,38) = 2.1;
% x(1,45:55) = ones(1,11)*3.8;
% N = length(x);%��ȡ��������
% t = (0:N-1);%��ʾʵ��ʱ��
% hold on
% plot(t,x,'g');%����ʱ����
% gausFilter = fspecial('gaussian',[1 3],3);
% v1 = ones(1,63)*gausFilter(1);
% L1 = diag(v1,1);
% v2 = ones(1,64)*gausFilter(2);
% L2 = diag(v2);
% v3 = ones(1,63)*gausFilter(3);
% L3 = diag(v3,-1);
% L = L1 + L2 + L3;
% u = rand(1,64)*0.15;
% x2 = L*x'+ u';
% plot(t,x2,'r');
% A=[-L;L;eye(64);-eye(64)];
% b=[-x2-0.15;x2-0.15;zeros(64,1);-ones(64,1)*12];
% x0=zeros(64,1);
% maxIter = 30;
% [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
%  plot(t,xkh','b');
%  
 
 
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