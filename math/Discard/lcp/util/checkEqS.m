function y=checkEqS(s1,s2)
I1=(s1>0);
I2=(s2>0);
y=any(xor(I1,I2));