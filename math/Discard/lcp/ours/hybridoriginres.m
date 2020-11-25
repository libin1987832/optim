function [x,err]=hybridoriginres(x0,nmax,M,q)
addpath('./predict');
% tol so that machine error
tol=1e-10;
% the number of the split method iteration
nf=5;
% error for x0
err=test_valid(M,q,x0);
while err>tol && index< nmax
    xkA = splitS(M,q,1,x0,nf);
    I=(xkA(:,nf-2)>0);
    % predict forumla eigen
    t=predict2(xkA(:,nf-2),xkA(:,nf-1),xkA(:,nf));
    t(t>0)=1;
    t(t<0)=0;
    cs=sum(I);
    cs1=sum(t);
    % compare the sum of the number(paper ask set is equivalent) 
    if cs1==cs
        % subspace step
        xs=subspacesearch(xkA(:,nf-2),M,q);
        % check for optimality
        err=test_valid(M,q,xs);
        % if err ok return function
        x=xs;
    else
        % if predit no exposed point ,contiune Guass-seidel
        x0=xkA(:,nf);
        x=x0;
        err=test_valid(M,q,x0);
    end
end
csA