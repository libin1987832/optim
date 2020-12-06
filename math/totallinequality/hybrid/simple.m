function [xkA, rpk] = simple(A, b, x0, n, rpk, nf, tol, options, type)
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
            maxit = 8;
            tou = 0.6;
            z = -rpk;
            z(z<0) = 0;
            bk = b + z;
            [x0,~] = IPGnnls(A, bk, x0, rpk, tol, tou,maxit);
            rpk = b - A * x0;
            xkA(:,i) = x0;
        end        
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