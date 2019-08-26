function [y,AI,bI]=Areax(x0,x1,A,b)
r0=b-A*x0;
I0=(r0>=0);

r1=b-x1;
I1=(r1>=0);

y=isequal(I0,I1);
AI=[];
bI=[];
if y==0
    AI=A(I0,:);
    bI=b(I0);
end