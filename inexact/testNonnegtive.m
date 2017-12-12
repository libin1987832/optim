addpath(genpath(pwd));
A=[2,1;-3,-1];
b=[5;-3];
A=-1*A;
b=-1*b;
%  [x,f0]=alg1(A,b,[0;100],0.0001)
% [x,fk]=fsearchx(A,b,[10;1],0.000001,0.0001)
 [x,fk]=inexact(b,A,[10;1],0.5,0.6)
% [mm,vv]=minValue(x,0.1,100,A,b,@fq,@constrain);
% mm
% ezplot('2*x+y-5',[0,10])
% hold on 
% ezplot('3*x+y-3',[0,10])