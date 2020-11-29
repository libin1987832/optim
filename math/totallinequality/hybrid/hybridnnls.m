function [xk,resvec,arvec,face1vec,face2vec,tf] = hybridnnls(A,b,x0,alpha,nf,maxit,options)
t=clock;

% stop criterion
tol = 1e-13;
[m,n] = size(A);
lsqrTol = 1e-13;
maxIter = maxit;
[rpk, normr, xmin, Ar, normKKT,face1,face2] = kktResidual(A, b, x0,[],1);
iter = 0;
% the residual vector
resvec = zeros(1,(maxit + 1)*(nf+1));
% the normal gradient 
arvec = zeros(1,(maxit + 1)*(nf+1));
% subspace minization
itersm = zeros(1,maxit + 1);
% face b-Ax
face1vec = zeros(1,(maxit + 1)*(nf+1));
% face x
face2vec = zeros(1,(maxit + 1)*(nf+1));
resvec(1) = xmin;
arvec(1) = normKKT;
face1vec(1) =face1;
face2vec(1) =face2;
indexsm = 0;
% flag 0-4 return lsqr flag
flag = 5;
xfA = zeros(n,nf);
xfA(:,1) = x0;
while normKKT > tol 
    iter = iter + 1;
    for i = 2:nf+1
    % AA indice of active set
    AA = (rpk>tol);
    % subspace
    AI = A(AA,:);
    bI = rpk(AA);
    [u,flag,relres,~,~,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA,2);
    xls = x0 + u;
    %[x,rpk] = projectfixed(A,b,xls,rpk,alpha);
    z = -rpk;
    z(z<0) = 0;
    bk = b + z;
    [x,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),xls,options);
    x0 = x;
    rpk = b - A * x;
     [rpk, normr, xmin, Ar, normKKT, face1, face2] = kktResidual(A, b, x ,rpk,1);
     resvec((iter-1)*(nf+1) + i-1) = normr;
    % record the value of the gradient function
    arvec((iter-1)*(nf+1) + i-1) = normKKT;
    face1vec((iter-1)*(nf+1) + i-1) = face1;
    face2vec((iter-1)*(nf+1) + i-1) = face2;    
    xfA(:,i-1) = x;
    end
    isSub = true;
   % isSub = strategy(A,b,x0,[],'PHA',iter,nf,rpk,xfA);
    %isSub = strategy(A,b,steplengthOrk,2,iter,nf,rpk,xfA);
    if isSub
        AA = (rpk>tol);
        RR = (x0>0);
        % subspace
        AI = A(AA,RR);
        bI = rpk(AA);
        u = lsqminnorm(AI,bI);
        p = zeros(n,1);
        p(RR) = u;
        [alpha, minf, knot] = arraySpiecewise(A,b,x0,p);
%         
% %         knoty = arrayfun(@(alpha) funmin(A,b,x0,u,alpha), knot);
% 
%         xa = [0:0.001:0.1];
% ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
%  pxy={};
% % pxy(1).X = knot;
% % pxy(1).Y = knoty;
% pxy(1).X = xa;
% pxy(1).Y = ya;
% figure
% hold on
% p1 = arrayfun(@(a) plot(a.X,a.Y),pxy);
% p1(1).Marker = 'o';
% %p1(2).Marker = '+';
% hold off
% r = b - A * x0;
% r(r<0)=0;
% Ar = A'*r;
% [alpha,Ar'*p]

        x0 = x0 + alpha*p;
        [rpk1, normr1, xmin1, Ar1, normKKT1, face12, face22] = kktResidual(A, b, xfA(:, end) , [], 1);
        [rpk, normr, xmin, Ar, normKKT, face1, face2] = kktResidual(A, b, x0 , [], 1);
            
        %u = AI\bI;
        %[u,flag,relres,~,resvec,lsvec,out] = lsqrm(AI,bI,lsqrTol,maxIter,[],[],zeros(n,1),A,b,x0,AA,3);
        indexsm = indexsm + 1;
%         xls = x0 + u;
%         x0 = xls;
    else
        x0 = xfA(:, end);
    end
    [rpk, normr, xmin, Ar, normKKT, face1, face2] = kktResidual(A, b, x0 , [], 1);
    % record the value of objection function
    resvec((iter-1)*(nf+1) + nf+1) = normr;
    % record the value of the gradient function
    arvec((iter-1)*(nf+1) + nf+1) = normKKT;
    face1vec((iter-1)*(nf+1) + nf + 1) = face1;
    face2vec((iter-1)*(nf+1) + nf + 1) = face2;
%    if iter > maxit || flag == 0
if iter > maxit
        break;
    end
    if flag < 5 || flag > 6 
%         if flag < 12
%             disp(['flag:',num2str(flag)]);
%         end
        flag = 5;
    end
end
xk = x0;
resvec = resvec(resvec>0);
% the normal gradient 
arvec = arvec(arvec>0);
% subspace minization
itersm = itersm(itersm>0);
% face b-Ax
face1vec = face1vec(face1vec>0);
% face x
face2vec = face2vec(face2vec>0);

tf = etime(clock,t);