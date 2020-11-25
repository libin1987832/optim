
% A=[1,1;-1,-1;-1,0;-6,-3];
% b=[1;1;0.5;2];

A=[1,1;-1,-1;-1,0];
b=[1;1;0.5];

[q,r]=qr(A);

xx1=[];
yy1=[];
xx2=[];
yy2=[];
nIter=3;
count=0;
z1=[];
z2=[];
z3=[];
for x=-2.5:0.1:0
    for y=-1.5:0.2:3
% x=1;
% y=-1;
       count=count+11;
        x0=[x;y];
        z=A*x0-b;
        z(z<0)=0;
        z1=[z1,z(1)];
        z2=[z2,z(2)];
        z3=[z3,z(3)];
         % in the expose face 
  %      x0=[-2.5;1.7]
        % not in
     end
 end
% express neton's method good
plot3(z1,z2,z3,'p','MarkerSize',10)
x1=xlabel('z1Öá');        
x2=ylabel('z2Öá');      
x3=zlabel('z3Öá');        
hold on 
% express FM method good
%plot(xx2,yy2,'bo');


