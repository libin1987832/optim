% the derive for (b-Ax)+
function detF=dFM(A,b,x1)
r=b-A*x1;
r(r<0)=0;
detF=A'*r;