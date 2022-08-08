n=100;
C=unidrnd(10,n,n);
[C,Y]=qr(C,0);
C;
dlmwrite('C1.txt', C, 'precision', '%5f', 'delimiter', '\t')