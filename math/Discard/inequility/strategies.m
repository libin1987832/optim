% parameter: A b x0 n obvious
% Qn is the decomposition for A if Dax method
% Qn may be steplength if the gradient method
% rpk is b-Ax
function isSub = strategies(A,b,Qn,type,iter,nf,rkp,xA)
[m,n] = size(A);
isSub = false;
eIter = 2;
con1 = 0.95;
con2 = 0.2;
diff = 100*eps;
tol = 1e-13;
Daxiter = floor(max(33,(m+n)/4)/nf);
switch upper(type)
    case 'GHA'
        rk=rkp;
        AA=find(rk>1e-15);
        qrkn=A(AA,:)'*rk(AA);
        qrkn=A*qrkn;
        rkn=rkp-eIter*Qn*qrkn;
        ssign=sum(~xor(rk>tol, rkn>tol));
        if ssign==m 
            isSub = true;
        end
    case 'DHA'
        if mod(iter,Daxiter) == 0
            isSub = true;
        end
    case 'PHA'
        rk=rkp;
        AA=find(rk>1e-15);
        qrkn=Qn(AA,:)'*rk(AA);
        qrkn=Qn*qrkn;
        rkn=rkp-eIter*qrkn;
        ssign=sum(~xor(rk>tol, rkn>tol));
        if ssign==m 
            isSub = true;
        end
    case 'CHA'
        x0 = xA(:,end - 1);
        xpn = xA(:,end - 2);
        xk = xA(:,end);
        p1v=x0-xpn;
        p1=p1v'*p1v;
        p2v=xk-x0;
        p2=p2v'*p2v;
        roup=p2/p1;
        if roup<con1 && roup>con2
            isSub = true;
        end
     case 'RHA'
        x0 = xA(:,end-1);
        xk = xA(:,end);
        rpk0 = b - A * x0;
        rpk = b - A * xk;
        r0 = rpk0;
        rk = rpk;
        r0(r0<0) = 0;
        rk(rk<0) = 0;
        %|| L(xk,zk) ||
        L1 = r0'*r0;
        %||L(xk+1,zk+1)||^2
        L2 = rk'*rk;
        z0=-rpk0;
        z0(z0<0)=0;
        % -(b-Axk+1)-zk
        LMv=-rpk-z0;
        LM=LMv'*LMv;
        %L(xk,zk)-L(xk+1,zk)
        L1m=L1-LM;
        %L(xk+1,zk)-L(xk+1,zk+1)
        L2m=LM-L2;
        %[L(xk,zk)-L(xk+1,zk)]-[L(xk+1,zk)-L(xk+1,zk+1)]
        LL=L1m-L2m;
        
        % active set r0
        AAz=(rkp>tol);
        %Au=-(b-Ax_{k+1})+(b-Axk)
        Au=-rpk+rpk0;
        % Au(AAz)
        AuAA=Au(AAz);
        %|| Au(AAz) ||^2
        LLA=AuAA'*AuAA;
        if abs(LLA-LL) < diff
            isSub = true;
        end
    otherwise
        error(message('MATLAB:InvalidOp'))
end
