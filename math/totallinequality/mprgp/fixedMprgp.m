function [x0, resvec, arvec, face1vec, face2vec, tf] = fixedMprgp(A,b,x0,param)
t=clock;
[m,n]=size(A);
index = 1;
% the residual vector
resvec = zeros(1, param.fixed_maxit + 1);
% the normal gradient 
arvec = zeros(1, param.fixed_maxit + 1);
face1vec = zeros(1,(param.fixed_maxit + 1));
% face x
face2vec = zeros(1,(param.fixed_maxit + 1));
[r, normr, xmin, Ar, KKT, face1, face2] = kktResidual(A, b, x0 , [],1);
resvec(index) = normr;
arvec(index) = KKT;
face1vec(index) = face1;
face2vec(index) = face2;

  
while KKT > param.fixed_tol && index < param.fixed_maxit
    x1=MPRGP(A, b, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter);
    x0=x1;
    [r, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x0 , [], 1);
    fprintf('iter=%d, normr=%g, Ar=%g\n', index, normr, KKT);
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
