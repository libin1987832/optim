n=100;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
C;
dlmwrite('C2.txt', C, 'precision', '%5f', 'delimiter', '\t')

%% ���ɶԳ���������A 
a=100;
b=-50;
z=zeros(n,1);
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
z
Q=diag(z);% �����䣨a��b���ڲ������ȷֲ���nά����
A=C'*Q*C; %��ʼ����A
dlmwrite('A2.txt', A, 'precision', '%5f', 'delimiter', '\t')
%% ���ɶԳ���������B
C=unidrnd(10,n,n);
I=eye(n);
B=C'*C+I; %��ʼ����B
dlmwrite('B2.txt',B, 'precision', '%5f', 'delimiter', '\t')