function [x,Out]=dARGauss_Seidel(A,b,opts)

% adaptive randomized  Gauss-Seidel method for solving linear systems
%              Ax=b
% 
% the  initial x=0, hence A'r_0=A'(Ax_0-b)=-A'b
%
%Input: the coefficent matrix A, the vector b and opts
%opts.p: the choice for setting p
%opts.TOL: the stopping rule
% 
%.....
%
%Output: the approximate solution x and Out
% Out.error: the relative iterative residual \|Ax_k-b\|/\|Ax_0-b\|
% Out.iter: the total number of iteration
% ....
%
% Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn
%


[m,n]=size(A);

%% setting some parameter
flag=exist('opts');

if (flag && isfield(opts,'Max_iter'))
   Max_iter=opts.Max_iter;
else
    Max_iter=200000;
end

if (flag && isfield(opts,'TOL'))
   TOL=opts.TOL;
else
    TOL=10^-6;
end

if (flag && isfield(opts,'strategy'))
   strategy=opts.strategy;
else
    strategy=1;
end

if (flag && isfield(opts,'xstar'))
    xstar=opts.xstar;
    normxstar=norm(xstar);
   relative_error=1;
else
    relative_error=0;
end

if (flag && isfield(opts,'initialx'))
   initialx=opts.initial;
else
    initialx=zeros(n,1);
end

%%
x=initialx;


colunmnormA=sum(A.^2,1);
      pro1=colunmnormA/sum(colunmnormA);
     residualvector=-b;
     cumsumpro=cumsum(pro1)';
 if strategy==1
    l1=sum(repmat(cumsumpro,1,Max_iter)<repmat(rand(1,Max_iter),n,1),1)+1; 
 else
     B=A'*A;
 end


%%
sNresidual_r=norm(b)^2;
At_r=A'*b;
snormb=norm(b)^2;
Axberror=[];

%% 
stopc=0;
iter=0;
iter1=0;
while ~stopc
    iter=iter+1;
    iter1=iter1+1;
%     if iter1>m
%        iter1=1;
%        l1=sum(cumsumpro<rand(m,1),2)+1; %一次选取多个指标
%     end
    
    %% setting the probabilities for choosing the index
    switch strategy
        case 2
            p=opts.p;
            pnormAx_b=power(abs(At_r)./sqrt(colunmnormA)',p);
            prob=(pnormAx_b/sum(pnormAx_b))';
            cumsumpro=cumsum(prob);
    end
    %% selecting an index with (adaptive) probability
    if strategy==3
         [~,l]=max(abs(At_r./sqrt(colunmnormA)')); % choosing the most coherence one;
    else if strategy==2
         %l=randsrc(1,1,[1:n;prob]);
         l=sum(cumsumpro<rand)+1;
        end
    end
   
    %% update x
    %%%%%%%%%%%%%% the classical Gauss-Seidel method
    if strategy==1
        l=l1(iter1);
        midv=A(:,l)'*residualvector/colunmnormA(l);
         x(l)=x(l)-midv;
        residualvector=residualvector-midv*A(:,l);
%         if relative_error
%            %  error1=norm(x-xstar)/normxstar;
%            %  Axberror=[Axberror,error1];
%              if error1<TOL | iter>Max_iter
%                  stopc=1;
%              end
%         else
       % normres=norm(residualvector);
       % Axberror=[Axberror,normres/sqrt(snormb)];
       % if (normres^2/snormb)<TOL | iter>Max_iter-1
       if iter>Max_iter-1
        stopc=1;
       % end
        end
    else
         %%%%%%%%%%%%%% the adaptive Gauss-Seidel method
    x(l)=x(l)-(At_r(l))/colunmnormA(l); 
    %% update \|r_k\|^2
    sNresidual_r=sNresidual_r-((At_r(l))^2/colunmnormA(l));
    
    %% storing the residual error
   % Axberror=[Axberror,sqrt(sNresidual_r)/sqrt(snormb)];
   
    %% checking whether the stopping rules is satisfied
  %  if (sNresidual_r/snormb)<TOL | iter>Max_iter
  if iter>Max_iter
        stopc=1;
    else
          %% update A'r_k=A'(Ax_k-b)
          %Al=A(:,l);
          %AtAl=A'*Al;
    At_r=At_r-((At_r(l))/colunmnormA(l))*B(:,l);
    end
    end
    
end
%% setting Output
Out.error=Axberror;
Out.iter=iter;

end

