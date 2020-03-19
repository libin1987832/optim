% paper rand [1,1] example
%<<<<<<< HEAD
m=300;
ratio=0.6;
n=ceil(ratio*m);
 A=2*rand(m,n)-1;
 b=2*rand(m,1)-1;
% load('ddf.mat')
save('ddf.mat','A','b')
x0=zeros(n,1);

[xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b);
[xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);
% [xkArr1m,xkArr1n]=size(xkArr1);
% [xkArr2m,xkArr2n]=size(xkArr2);
% beginRow=18;
% plx1=beginRow:xkArr1m;
% ply1=xkArr1(plx1,n+1)';
% plx2=beginRow:xkArr2m;
% ply2=xkArr2(plx2,n+1)';
% plot(plx1(xkArr1(plx1,n+2)'==1),ply1(xkArr1(plx1,n+2)'==1),'ro');
%  hold on 
% plot(plx1(xkArr1(plx1,n+2)'==0),ply1(xkArr1(plx1,n+2)'==0),'r.');
% plot(plx2(xkArr2(plx2,n+2)'==1),ply2(xkArr2(plx2,n+2)'==1),'bo');
% plot(plx2(xkArr2(plx2,n+2)'==0),ply2(xkArr2(plx2,n+2)'==0),'b-');

% for m=1000:1000:5000
%     for ratio=0.1:0.3:0.7
% % m=300;
% % ratio=0.5;
%=======
% m=300;
% ratio=0.2;
%>>>>>>> a73cbedf22b6ac1a723380c9f9c07baa55399595
% n=ceil(ratio*m);
% A=2*rand(m,n)-1;
% b=2*rand(m,1)-1;
% x0=zeros(n,1);
% 
% [xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b);
% [xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);


% for m=1000:2000:5000
%     for ratio=0.1:0.3:0.7
% % m=300;
% % ratio=0.5;
% n=ceil(ratio*m);
% A=2*rand(m,n)-1;
% b=2*rand(m,1)-1;
% x0=zeros(n,1);
% 
% [xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b);
% [xk2,fk2,xkArr2,countF2,countN2]=hybrid2(x0,A,b);
%     end
% end
% [xk2,fk2,xkArr2,countF2,countN2]=hybrid4(x0,A,b);

