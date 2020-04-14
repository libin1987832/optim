% take the positive (b-Ax) 
function [y,AI,bI]=Areax(x0,x1,A,b)
% take the positive of A
r0=b-A*x0;
I0=(r0>=0);
% take the positive of b
r1=b-x1;
I1=(r1>=0);

% check if equal and return
y=isequal(I0,I1);
AI=[];
bI=[];
if y==0
    AI=A(I0,:);
    bI=b(I0);
end