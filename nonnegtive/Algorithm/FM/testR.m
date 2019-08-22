b=[1,1,2]';
A=[2,0;-2,0;0,1];
[Q,R]=qr(A);
x0=[-10;-10];
display=[];
for i=0:22
    [xk,y1,y2,fk,fm,fr]=FM(x0,Q,R,A,b);
    ro=(fr-fm)/(fr-fk);
    display=[display [ro;x0;y1]];
    x0=xk;
end
display