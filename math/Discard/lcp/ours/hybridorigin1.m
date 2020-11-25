% index: statistic split iterations indexN: statistic subspace minimzaiton
% allValue:record all value,[split predict xs(zeors for no exposed point)]
function [x,errA,index,indexN]=hybridorigin1(x0,nmax,nf,M,q,iDLU,iDL)
tol=1e-10;
err=test_valid(M,q,x0);
index=0;
indexN=0;
count=nmax*nf;
errA=cell(6,nmax);
tol_rel  = 1e-5;
tol_abs  = 1e-10;
while err>tol && index< count
    index=index+nf;
%     xkA= splitS_our(M,q,1,x0,nf);
%     xkA= splitS_ourres(M,q,1,x0,nf);
    if index/nf<3
    [ynf0 err iter flag convergence msg] =  pgs(M,q, x0, nf-2, tol_rel, tol_abs, true);
    [ynf1 err iter flag convergence msg] =  pgs(M,q, ynf0, 1, tol_rel, tol_abs, true);
    [ynf2 err iter flag convergence msg] =  pgs(M,q, ynf1, 1, tol_rel, tol_abs, true);
    else
    ynf0=Gauss_Siedel(iDLU,iDL,q,x0,nf-2);
    ynf1=Gauss_Siedel(iDLU,iDL,q,x0,1);
    ynf2=Gauss_Siedel(iDLU,iDL,q,x0,1);
    end
    xkA=[ynf0 ynf1 ynf2];
%     I=(xkA(:,3)>0);
%     I=(xkA(:,nf-2)>0);
%     t=predict2(xkA(:,nf-2),xkA(:,nf-1),xkA(:,nf));
    t=predict2(xkA(:,1),xkA(:,2),xkA(:,3));
%     t(t>0)=1;
%     t(t<0)=0;
%     cs=sum(I);
%     cs1=sum(t);
      [cs1,cs]=checkEq(t,xkA(:,3));
    errA(1,index/nf)={cs};
    errA(4,index/nf)={xkA(:,3)}; 
    if cs1==cs
        % subspace step
%         xs=subspacesearch(xkA(:,nf-2),M,q);
%         allValue=[allValue xs];
        % check for optimality
        xs=subspacesearch(xkA(:,3),M,q);
        err=test_valid(M,q,xs);
        errA(3,index/nf)={'subspace'}; 
        errA(5,index/nf)={xs}; 
        indexN=indexN+1;
    else
        err=test_valid(M,q,xkA(:,3));
        errA(3,index/nf)={'guass-seidel'}; 
        errA(5,index/nf)={xkA(:,3)}; 
    end
    errA(2,index/nf)={err};
%      x0=xkA(:,nf);
    x0=xkA(:,3);
    errd=test_valid(M,q,xkA(:,3));
    errA(6,index/nf)={errd}; 
%      disp(['hy err:',num2str(err)]);
end
if index< count && errd>err
x=xs;
else
% x=xkA(:,nf);
x=xkA(:,3);
end
err=test_valid(M,q,x);