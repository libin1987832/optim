% A=[1 0 0;2 1 0;3 0 1];
% B=[1 0 0;0 1 0;0 4 1];
% inv(A)*inv(B)
% A*B
% B*A 
A=[1;-1];
b=[-1;2];
[Q,R]=qr(A);
x0=1;
r=b-A*x0;
[xk1,rk1]=FixedM(x0,Q,R,A,b,r)
[xk2,rk2]=FixedM(xk1,Q,R,A,b,rk1)
[xk3,rk3]=FixedM(xk2,Q,R,A,b,rk2)
[xk4,rk4]=FixedM(xk3,Q,R,A,b,rk3)
r(r<0)=0;
rk1(rk1<0)=0;
rk2(rk2<0)=0;
rk3(rk3<0)=0;
rk4(rk4<0)=0;
rka=[r,rk1,rk2,rk3,rk4];
plot(rka(1,:),rka(2,:),'o');
hold on
line([0,1],[0,1],'LineStyle','--');
line([1,-1],[-1,1],'LineStyle','--');
line([-1,3],[0,0],'LineStyle','-');
line([0,0],[-1,3],'LineStyle','-');
axis([-1 3 -1 3])
