function [x0, resvec, arvec, face1vec, face2vec, tf] = fixedMprgp(A,b,x0,param)
t=clock;
[m,n]=size(A);
index = 1;
% the residual vector
resvec = zeros(5, param.fixed_maxit + 1);
% the normal gradient 
arvec = zeros(1, param.fixed_maxit + 1);
face1vec = zeros(1,(param.fixed_maxit + 1));
% face x
face2vec = zeros(1,(param.fixed_maxit + 1));
[rpk, normr, xmin, Ar, KKT, face1, face2] = kktResidual(A, b, x0 , [],1);
resvec(index) = normr;
arvec(index) = KKT;
face1vec(index) = face1;
face2vec(index) = face2;

options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options = optimoptions('Algorithm','interior-point','TolX',1e-13)
options.Display = 'off';
% options.StepTolerance = 1e-13;
options.OptimalityTolerance = 1e-15;
% options.ConstraintTolerance = 1e-13;
options.MaxIterations = 600;
 t0=1;
 z0=0;
 ATA=A'*A;

 debug=0;% 1:MPRGQ 2:MPRG 4:LS
 method1 = 1;
 method2 = 1;
while KKT > param.fixed_tol && index < param.fixed_maxit
    % update z    
    z = A*x0-b;
    z(z<0) = 0;
%     bk=b+z;
    % update t1
    t1=(1+sqrt(1+4*t0^2))/2;
    zn = z+(t0-1)/(t1)*(z-z0);
   % z0=z;
    z0=zn;
    t0=t1;
    if method2 == 1
     bk = b + zn;
    end
    if method2 == 2
      bk = b +z;
    end
     if method1 == 1  
        [x1,rpk]=MPRGPQ(ATA, A'*bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter,bk'*bk);
     end
     if method1 == 3
      [x1,rpk]=MPRGP(A, bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter);
     end
    if method1 == 4
    %内置使用MATLAB自带的算法
     [x1,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
    end
    x0=x1; 
    if debug
    [rpk, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x1 , [], 1);
    fprintf('iter=%d, normr=%g, KKT=%g,xminkr=%g\n', index, normr, KKT,xmin);     
    end
    index=index+1;
    resvec(index) = normr;
    arvec(index) = KKT;

    face1vec(index) = face1;
    face2vec(index) = face2;
end
resvec = resvec(resvec>0);
% the normal gradient 
arvec = arvec(arvec>0);
% face b-Ax
face1vec = face1vec(face1vec>0);
% face x
face2vec = face2vec(face2vec>0);

tf = etime(clock,t);
% fprintf('exact run time cost:%f\n',etime(clock,t));

%        if index == 15
%             F = x0 > param.mprgp_Ftol;
%             AF=ATA(F,F);
%             L = ichol(AF,struct('michol','on'));
%             [x2,fl2,rr2,it2,rv2] = pcgself(AF,A(:,F)'*bk,1e-2,3,L,L',x0(F));
%             x00=x0;
%             x00(F)=x2;
%             [rpk, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x00 , [], 1);
%             fprintf('iter=%d, normr=%g, KKT=%g,xminkr=%g\n', index, normr, KKT,xmin); 
%             [x1,rpk]=MPRGPQ(ATA, A'*bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter,bk'*bk);
%              norm(x2-x1(F))
%         end
