n=100;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
C;
dlmwrite('C2.txt', C, 'precision', '%5f', 'delimiter', '\t')

%% 生成对称正定矩阵A 
a=100;
b=-50;
z=zeros(n,1);
for i=1:n
    z(i)=b+((i-1)*(a-b))/(n-1);
end
z
Q=diag(z);% 在区间（a，b）内产生均匀分布的n维向量
A=C'*Q*C; %初始矩阵A
dlmwrite('A2.txt', A, 'precision', '%5f', 'delimiter', '\t')
%% 生成对称正定矩阵B
C=unidrnd(10,n,n);
I=eye(n);
B=C'*C+I; %初始矩阵B
dlmwrite('B2.txt',B, 'precision', '%5f', 'delimiter', '\t')