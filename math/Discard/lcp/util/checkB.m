% is in the feasible
function [b,s]=checkB(M,q,x0)
xn=any(x0<0);% if negative
f0=M*x0+q;
fn=any(f0<0); % if negative
 su=sum(x0.*f0); % the complement condition
if (xn+fn)==0 && su <1e-5
    s=1;
    A=(x0<1e-5);
    C=(f0>1e-5);
    bf=A & C;
    b=find(bf);
else
    s=0;
    b=[];
end