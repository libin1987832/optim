% A=[1,1;-1,-1;-1,0;-6,-3];
% b=[1;1;0.5;2];
% x0=[-1/2;1/3];

qI1=[1;-1;-1;-6];
q1=qI1/sqrt(39);
qI2=[19;-19;20;-3];
q2=qI2/sqrt(1131);
Q=[q1,q2];
Q1=qI1*qI1'*1131;
Q2=qI2*qI2'*39;
rT=diag([39*1131,39*1131,39*1131,39*1131])-(Q1+Q2)*diag([1,1,0,0]);
rts=sym(rT/(39*1131));
for i=-2:0.01:-1/2
    for j=-1-i:0.1:1-i
        x0=[i;j];
        if 2-[-6,-3]*x0>0
            continue;
        end
        x0s=sym(x0);
        % rtsb=rts*[1;1;0.5;2]
        % rtsA=rts*[1,1;-1,-1;-1,0;-6,-3]
        % rtsb-rtsA*x0s
        rs=[1;1;0.5;2]-[1,1;-1,-1;-1,0;-6,-3]*x0s;
        fm=rts*rs;
        fmi=fm.*rs;
        imin=min(fmi);
        if imin <0
            x0s
            [rs,fm]
        end
    end
end

% fm.*rs;
% r=[3/2;1/2;-1/4;-7/4];
% rs=sym(r);
% rst=v*j^10*inv(v)*rs
% func(47/29,18/29,1/2,-1/2-1/100)
function NK=NMatrix(A,b,x)
r=b-A*x;
o=ones(size(r));
o(r<0)=0;
o(r>0)=1;
NK=diag(o);
end
function x2=func(a1,a2,b,x1)
x2=-(b+a1*x1)/a2;
end