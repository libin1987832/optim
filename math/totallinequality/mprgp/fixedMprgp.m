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
 t0=1;
 z0=0;
 ATA=A'*A;

 debug=1;% 1:MPRGQ 2:MPRG 3:LS
 if debug ==2
%      L=spdiags(sqrt(1./(spdiags(ATA,0))),0,n,n);
    L = ichol(ATA)';
    H=inv(L)';
    L=H';
 %  L=diag(sqrt(1./diag(ATA)));
     LA=A*H;
      ATA=H*ATA*H';
      param.mprgp_a = 1/norm(ATA, inf)*30;
      y0=L\x0;
 end
while KKT > param.fixed_tol && index < param.fixed_maxit
    % update z    
    z = A*x0-b;
    z(z<0) = 0;
    bk=b+z;
    % update t1
    t1=(1+sqrt(1+4*t0^2))/2;
    zn = z+(t0-1)/(t1)*(z-z0);
   % z0=z;
    z0=zn;
    t0=t1;
     bk = b + zn;
     bk = b +z;
     if debug == 1
    [x1,rpk]=MPRGPQ(ATA, A'*bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter,bk'*bk);
     end
    if debug == 2
    [y1,rpk]=MPRGPQ(ATA, LA'*bk, y0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter,bk'*bk);
   % y1'*ATA*y1-2*(LA*y1)'*bk+bk'*bk
   % [yy1,fl0,rr0,it0,rv0]=pcg(ATA,LA'*bk,1e-10,3); rv0=norm(LA'*bk-ATA*yy1,2) yy1'*ATA*yy1-2*(LA*yy1)'*bk+bk'*bk
    x1=L*y1;
    y0=y1;
     end
     if debug == 3
      [x1,rpk]=MPRGP(A, bk, x0, param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, param.mprgp_maxIter);
     end
    if debug == 4
    %内置使用MATLAB自带的算法
     [x1,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
    end
     x0=x1;
%     rpk = b-A*x1;
%     z=-rpk;
%     z(z<0)=0;
%     bk = b + z;
%     z0=z;
%      [x1,f1,residual,exitflag,output,ff]=lsqlin(A,bk,[],[],[],[],zeros(n,1),Inf*ones(n,1),x0,options);
 
%     F = x0 > param.mprgp_Ftol;
%     xkkr = x0;
%     u = krylovkm(A(:,F),b,rpk,10);
%     xkkr(F)=xkkr(F)+u;
%     [rpkkr, normrkr, xminkr, Ar, KKTkr, face1, face2] = kktResidual(A, b, xkkr , [],1);
   
  

     
%     
         [rpk, normr, xmin,Ar, KKT,face1,face2] = kktResidual(A, b, x1 , [], 1);
         fprintf('iter=%d, normr=%g, KKT=%g,xminkr=%g\n', index, normr, KKT,xmin);     
  
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
