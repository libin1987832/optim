% FM+exposed change ����uIter FM ��û�иı������ ţ�ٷ�
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
uIter=33;
face=zeros(1,uIter);
xkArr=[];

countFM=0;
countNW=0;

beginNW=0;

if Ar<delt*rn || rn<delt
    rk=r;
    xk=x0;
    disp('input x is satisfied all constrain!') %ceases execution
end
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
        ArN=norm(A'*rk);
        rnN=norm(rk);
        if ArN<delt*rnN || rnN<delt
            break
        end
    end
    
    % the first is successise unchange the second is feasible the third machine error
    if all(face(2:end)) && ~isempty(IkN) && ArN>delt*rn && rnN>delt
        countNW=countNW+1;
        if countNW ==1
            beginNW=countFM;
        end
        % take the binary search step length
        [xk,rk,fk2,f0,lambe]=ssqr(xk,A,b);
        xkArr=[xkArr;[xk',fk2,1]];
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    r=rk;
    x0=xk;
end
fk=0.5*rk'*rk;
tf=etime(clock,t);
vk=sum(sign(rk));
disp(['%hybrid2 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' uIter:',num2str(beginNW)]);
%disp(['$',num2str(m),'\times ',num2str(n),'$&FMEF&(',num2str(countFM),',',num2str(countNW),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
%disp(['well1033&our&',num2str(vk),'&',num2str(rn),'&',num2str(Ar),'&(',num2str(countFM),',',num2str(countNW),')&',num2str(beginNW)]);