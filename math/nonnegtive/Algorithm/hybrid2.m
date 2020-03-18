% FM+exposed change 连续uIter FM 都没有改变则采用 牛顿法 
function [xk,fk,xkArr,countFM,countNW,Q]=hybrid2(x0,A,b)
t=clock;
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

countFM=0;
countNW=0;

beginNW=0;
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    I=find(r>tol);
    %FM algorithm
    statFM=statFM+uIter;
    IkN=I;
    for ii=1:uIter
        countFM=countFM+1;
        Ik=IkN;
        [xk,r0,rk,fk1,fm,fr]=FM(x0,Q,R,A,b);
         fk=fk1;
        xkArr=[xkArr;[xk',fk1,0]];
        IkN=find(rk>tol);
        face(ii)=isequal(IkN,Ik);
        x0=xk;
    end
    
    if all(face(2:end)) && ~isempty(IkN)
        %     if isequal(I,Ik1)
        %newtonalgorithm
        y=1;
        while y
            statSS=statSS+1;
            countNW=countNW+1;
            if beginNW <1
                beginNW=countFM;
            end
            [xk2,fk2,y]=ssqr2(xk,A,b);
             fk=fk2;
            xkArr=[xkArr;[xk',fk1,1]];
            if abs(fk2 - fk1) < 1e-7 || fk2 < fk1 
                xk=xk2;
                rk2=b-A*xk2;
                rk2(rk2<0)=0;
                rk=rk2;
                break;
            end
        end
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    r=rk;
    x0=xk;
end
fk=0.5*rk'*rk;
tf=etime(clock,t);
disp(['hybrid2 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM),' cpu:',num2str(tf),' uIter:',num2str(beginNW)]);
disp(['$',num2str(m),'\times ',num2str(n),'$&FMEF&(',num2str(statFM),',',num2str(statSS),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
