function isSub = strategies(A,b,Qn,iter,type,rkp,xA)
[m,n] = size(A);
isSub = false;
eIter = 3;
con = 0.8;
diff = 10*eps;
Daxiter = floor(max(33,(m+n)/4));
switch upper(type)
    case 'DHA'
        if mod(iter,Daxiter) == 0
            isSub = true;
        end
    case 'PHA'
        rk=rkp;
        AA=find(rk>ee);
        qrkn=Qn(AA,:)'*rk(AA);
        qrkn=Qn*qrkn;
        rkn=rkp-eIter*qrkn;
        ssign=sum(~xor(rk>ee,rkn>ee));
        if ssign==m ||ssign>m*0.99
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
        if roup<con
            isSub = true;
        end
     case 'RHA'
        x0 = xA(1);
        xk = xA(2);
        rp0 = b - A * x0;
        rpk = b - A * xk;
        r0 = rp0;
        rk = rpk;
        r0(r0<0) = 0;
        %|| L(xk,zk) ||
        L1 = r0'*r0;
        %||L(xk+1,zk+1)||^2
        L2 = rk'*rk;
        z0=-rkp0;
        z0(z0<0)=0;
        % -(b-Axk+1)-zk
        LMv=-rkp-z0;
        LM=LMv'*LMv;
        %L(xk,zk)-L(xk+1,zk)
        L1m=L1-LM;
        %L(xk+1,zk)-L(xk+1,zk+1)
        L2m=LM-L2;
        %[L(xk,zk)-L(xk+1,zk)]-[L(xk+1,zk)-L(xk+1,zk+1)]
        LL=L1m-L2m;
        
        % active set r0
        AAz=(rkp>ee);
        %Au=-(b-Ax_{k+1})+(b-Axk)
        Au=-rkp+rkp0;
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
