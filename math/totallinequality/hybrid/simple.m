function [xkA, rpk] = simple(A, b, x0, n, rpk, nf, tol, options, type)
display = true;
if display
rpkOld = rpk;
end
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
    case 'lsqr'
        for i = 2:nf+1            
            z = -rpk;
            z(z<0) = 0;
            bk = b + z;
            
            AA = (rpk>-tol);
            % subspace
            AI = A(AA,:);
            bI = rpk(AA);
            lsqrTol= 1e-5;
            maxIter =10;
            [u,flag,relres,iter,resvec,lsvec,out] = lsqrmx(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA,3);
            if ~out && all(x0 + u > 0) 
                x0 = x0 + u;
                
%                 rpk = b - A * xs;
                [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, x0, [], 1);
                if display
                    r = rpkOld; r(r < 0) = 0; g0 = -A' * r; pg = g0'*u; [minx,loc]=min(x0);
                    fprintf('simple(lsqmx lsqr): normr(%g),pg(%g),min(%g,%d)\n',normr,pg,minx,loc);
                end
            else
%                 u = AI\bI;
                    u = lsqminnorm(AI,bI);
                [aa, knot, retcode] = arraySpiece(A,b,x0,u,1e-5,30);
%                 aa = spiecewise(A,b,u,x0);
                xs = x0 + aa * u;
                xs(xs<0) = 0;
%                 rpk = b - A * xs;
                x0=xs;
                [rpk, normr, minx, g, normKKT, face1, face2] = kktResidual(A, b, xs, [], 1);
                if display
                    r = rpkOld; r(r < 0) = 0; g0 = -A' * r; pg = g0'*u; [minx,loc]=min(x0);
                    fprintf('simple(project): norm(%g),pg(%g),min(%g,%d),alpha(%g)\n',normr,pg,minx,loc,aa);
                end
            end
             xkA(:,i) = x0;
        end
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