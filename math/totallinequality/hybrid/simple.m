function [xkA, rpk] = simple(A, b, x0, n, rpk, nf, options, type)
xkA = zeros(n,nf+1);
xkA(:, 1) = x0;
switch type
    case 'lsqlin'
        for i = 2:nf+1
            z = -rpk;
            z(z<0) = 0;
            bk = b + z;
            [x, ~, ~, ~, ~, ~]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
            x0 = x;
            xkA(:,i) = x0;
            rpk = b - A * x;
        end
    case 'IPG'
        for i = 2:nf+1
            z = -rpk;
            z(z<0) = 0;
            bk = b + z;
            Ax = b - rpk;
            ATb = A' * bk;
            
            ATAx = A' * Ax;
            Ar = ATAx - ATb;
            d = x0./ATAx;
            p = - d .* Ar;
            alphaAll = - x0./p;
            alphak = min(alphaAll(alphaAll>0));
            if isempty(alphak)
                alphak = 1;
            end
            talphak = tou * alphak;
            Ap = A * p;
            Apn = Ap' * Ap;
            alphaStar = - (p' * Ar) / Apn;
            alpha = min(talphak, alphaStar);
            x0 = x0 + alpha * p;
            
            xkA(:,i) = x0;
        end
        %     xa = [0:0.0001:talphak];
        %     ya = arrayfun(@(alpha) func(A,b,x0,p,alpha), xa);
        %     ya2 = arrayfun(@(alpha) normr + alpha * 0.25 * g' * p , xa);
        %     ya3 = arrayfun(@(alpha) func(A,b,x0,p,alpha), allalpha);
        %     plot(xa,ya,'+',xa,ya2,'o',allalpha,ya3,'x');
        %     hold on
        
        %[alpha, knots, retcode] = arraySpiece(A,b,x0,p);
        x0 = x0 + alpha*p;
        
    case 'lsqrlsqlin'
end


% for i = 2:nf+1
%     % AA indice of active set
%     AA = (rpk>tol);
%     % subspace
%     AI = A(AA,:);
%     bI = rpk(AA);
%     [u,flag,relres,~,~,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA,2);
%     xls = x0 + u;
%     %[x,rpk] = projectfixed(A,b,xls,rpk,alpha);
%     z = -rpk;
%     z(z<0) = 0;
%     bk = b + z;
%     [x,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),xls,options);
%     x0 = x;
%     rpk = b - A * x;
%     [rpk, normr, ~, g, normKKT, face1, face2] = kktResidual(A, b, x ,rpk,1);
%     resvec((iter-1)*(nf+1) + i-1) = normr;
%     % record the value of the gradient function
%     arvec((iter-1)*(nf+1) + i-1) = normKKT;
%     face1vec((iter-1)*(nf+1) + i-1) = face1;
%     face2vec((iter-1)*(nf+1) + i-1) = face2;
%     xfA(:,i-1) = x;
% end