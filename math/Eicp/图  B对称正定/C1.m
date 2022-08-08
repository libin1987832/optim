n=100;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
C;
dlmwrite('C2.txt', C, 'precision', '%5f', 'delimiter', '\t')
C=unidrnd(10,n,n);
I=eye(n);
B=C'*C+I; %≥ı ºæÿ’ÛB
dlmwrite('B3.txt',B, 'precision', '%5f', 'delimiter', '\t')