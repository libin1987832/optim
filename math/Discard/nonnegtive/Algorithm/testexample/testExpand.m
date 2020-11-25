A=[1,2,4,5;1,1,2,3;4,6,7,9;1,6,7,9];
B=A'*A
lamd=eig(B);
syms t
pt=(t-lamd(1))*(t-lamd(2))*(t-lamd(3))*(t-lamd(4))
pt1=expand(pt)
pt2=subs(pt1,t,B)
pt2*A'