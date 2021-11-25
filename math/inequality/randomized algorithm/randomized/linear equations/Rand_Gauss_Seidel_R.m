function [x,normr,index]= Rand_Gauss_Seidel_R(A, r, At_r,maxit,B,Acol,p,alpha)
%% Ax=b Ar probilities
%%
[m, n] = size(A);


%% º∆À„≤–≤Ó

x = zeros(n,1);

pnormAx_b=(((At_r).^2)./Acol);
% pnormAx_b=Acol;
prob=(pnormAx_b/sum(pnormAx_b));

cumsumpro=cumsum(prob);
normr=[norm(A'*r)];
index=[];
for i = 1:maxit
    pickedj=sum(cumsumpro<rand)+1;
%     [~,pickedj]=max(abs(At_r./sqrt(Acol)));
    index=[index,pickedj];
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
 %   inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;
  %  At_r =A'*r;
    At_r = At_r - inc*B(:,pickedj);
%      pnormAx_b=(((At_r).^2)./Acol);
% 
% %     pnormAx_b=(abs(At_r)./sAcol);
    pnormAx_b=power((abs(At_r)./sqrt(Acol)),p);
    prob=(pnormAx_b/sum(pnormAx_b));
    cumsumpro=cumsum(prob);
       normr=[normr,norm(A'*r)];
end

end
