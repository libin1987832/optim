% index: statistic split iterations indexN: statistic subspace minimzaiton
% allValue:record all value,[split predict xs(zeors for no exposed point)]
function [x,err,index,indexN]=PA(x0,nmax,nf,M,q)
addpath('./predict');
tol=1e-10;
err=test_valid(M,q,x0);
index=0;
indexN=0;
csA=[];
count=nmax*nf;
tol_rel  = 1e-5;
tol_abs  = 1e-10;

while err>tol && index< count
    index=index+nf;
    [ynf0 err iter flag convergence msg] =  pgs(M,q, x0, nf-2, tol_rel, tol_abs, true);
    [ynf1 err iter flag convergence msg] =  pgs(M,q, ynf0, 1, tol_rel, tol_abs, true);
    [ynf2 err iter flag convergence msg] =  pgs(M,q, ynf1, 1, tol_rel, tol_abs, true);
    [a,xpf,result]=projectedsearch(ynf2-x0,x0,M,q);
%      [a,xpf,result]=projectedsearch(xkA(:,nf)-x0,x0,M,q);
     I=(xpf>0);
%     MII=M(I,I);
%     qI=q(I);
%    t=predict(MII,xpf(I),qI,4);
%     t=predict2(xkA(:,nf-2),xkA(:,nf-1),xkA(:,nf)); %mathc splitS_our
 t=predict2(ynf0,ynf1,ynf2);
    t(t>0)=1;
    t(t<0)=0;
    cs=sum(I);
    cs1=sum(t);
    csA=[csA,cs1]; 
    if cs1==cs
        % subspace step
        xs=subspacesearch(xpf,M,q);
         [a,xk]=projectedsearch(xs-xpf,xpf,M,q);
        x0=xk;
        indexN=indexN+1;
    else
        x0=xpf;
    end
    % check for optimality
    err=test_valid(M,q,x0);
    disp(['pa err:',num2str(err)]);
end
x=x0;
csA;