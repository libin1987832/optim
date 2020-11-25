function [num1,num2]=checkEq(a1,a2)
t=(a1>0);
I=(a2>0);
num2=sum(I);
num1=sum(t);
        