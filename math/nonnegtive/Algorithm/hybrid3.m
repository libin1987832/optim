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
    
    if all(face(2:end)) && ~isempty(IkN)
        
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    r=rk;
    x0=xk;
end
fk=0.5*rk'*rk;
disp(['AT(b-A*x)+:',num2str(Ar),' fk:',fk,' ssqr:',num2str(statSS),' FM:',num2str(statFM)]);