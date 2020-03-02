% FM+exposed change 连续uIter FM 都没有改变则采用 若牛顿法非奇异则继续计算 否则直接不采用牛顿法
function [xk,fk,xkArr]=hybrid3(x0,A,b)
tol=0;
%compute hybrid uIter
[m,n]=size(A);
%FM need a qr decompose
[Q,R]=qr(A);
r=b-A*x0;
r(r<0)=0;
%condition for terminate
Ar=norm(A'*r);
rn=norm(r);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
statFM=0;
statSS=0;
uIter=3;
face=zeros(1,uIter);
xkArr=[];
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    I=find(r>tol);
    %FM algorithm
    statFM=statFM+uIter;
    IkN=I;
    for ii=1:uIter
        Ik=IkN;
        [xk,r0,rk,fk1,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk1,0]];
        IkN=find(rk>tol);
        face(ii)=isequal(IkN,Ik);
        x0=xk;
    end
    [xk,fk,y,isP]=ssqr3(x0,A,b);
    if isP
        if y
            disp('last step');
        else
            disp('FM may not reach the right face');
        end
    else
       %  disp("face is not positon");
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    r=rk;
    x0=xk;
end
fk=0.5*rk'*rk;
disp(['AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM)]);