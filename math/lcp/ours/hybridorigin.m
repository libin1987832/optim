% index: statistic split iterations indexN: statistic subspace minimzaiton
% allValue:record all value,[split predict xs(zeors for no exposed point)]
function [x,err,index,indexN,allValue]=hybridorigin(x0,nmax,nf,M,q)
addpath('./predict');
addpath('../util');
tol=1e-10;
err=test_valid(M,q,x0);
index=0;
indexN=0;
csA=[];
allValue=[];
while err>tol && index/nf< nmax
    index=index+nf;
    cmax=0;
    [xkA,cmax] = splitS(M,q,1,x0,nf,cmax);
    I=(xkA(:,nf-2)>0);
    t=predict2(xkA(:,nf-2),xkA(:,nf-1),xkA(:,nf));
   % allValue=[allValue xkA t];
    allValue=[allValue xkA];
    t(t>0)=1;
    t(t<0)=0;
    cs=sum(I);
    cs1=sum(t);
    csA=[csA,cs1]; 
    if cs1==cs
        % subspace step
        xs=subspacesearch(xkA(:,nf-2),M,q);
%         allValue=[allValue xs];
        % check for optimality
        err=test_valid(M,q,xs);
        %x0=xk;
        indexN=indexN+1;
        x=xs;
    else
%        x0=xpf;
        x=xkA(:,nf);
        err=test_valid(M,q,x0);
        allValue=[allValue zeros(size(x0))];
    end
     x0=xkA(:,nf);
end
csA