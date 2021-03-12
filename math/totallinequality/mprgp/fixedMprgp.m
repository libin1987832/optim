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
 t0=0;
 z0=0;
 ATA=A'*A;
 L=spdiag(sqrt(1./diag(ATA)));
  ATA=L*ATA*L;
while KKT > param.fixed_tol && index < param.fixed_maxit
    
        z = -rpk;
    z(z<0) = 0;
    zn = z+(t0-1)/(t0+1)*(z-z0);
    t0=(1+sqrt(1+4*t0^2))/2;
    z0=zn;
     bk = b + zn;
    [x1,rpk]=MPRGPQ(ATA, L*A'*bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter);
    x1=L/x1;
%     [x1,rpk]=MPRGP(A, bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter);
    [rpk, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x1 ,[], 1);
%    fprintf('iter=%d, normr=%g, KKT=%g,xminkr=%g\n', index, normr, KKT,xmin);     
    x0=x1;



%        z = -rpk;
%     z(z<0) = 0; 
%      bk = b + z;

    %内置使用MATLAB自带的算法
%     [x2,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
%       [rpk, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x2 , [], 1);
%       fprintf('iter=%d, normr=%g, KKT=%g,xminkr=%g\n', index, normr, KKT,xmin);     
%     F = x0 > param.mprgp_Ftol;
%     xkkr = x0;
%     u = krylovkm(A(:,F),b,rpk,10);
%     xkkr(F)=xkkr(F)+u;
%     [rpkkr, normrkr, xminkr, Ar, KKTkr, face1, face2] = kktResidual(A, b, xkkr , [],1);
   
  
    
      
    
 %   [r, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x0 , rpk, 1);
%      if KKT > KKTkr
%          x0=xkkr;
%      end
  %    KKTkr = KKT;
%   fprintf('iter=%d, normr=%g, KKT=%g,krylovKKt=%g,xminkr=%g\n', index, normr, KKT,KKTkr,xminkr);
    
  
  
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
