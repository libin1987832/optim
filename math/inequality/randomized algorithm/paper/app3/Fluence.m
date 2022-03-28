function [x]=Fluence(A,B,u,l,k,r,maxIter,nnls)
[m1,n]=size(A);
[m2,n]=size(B);
H=[1/sqrt(m2)*B;sqrt(r/m1)*A];
w=0;
iter=0;
while 1
    iter=iter+1;
dw=[1/sqrt(m2)*l;sqrt(r/m1)*(u+w)];
f = -H'*dw;
HTH=H'*H;
x0=zeros(n,1);
lb=zeros(n,1);
ub=inf*ones(n,1);
if strcmp(nnls,'quadprog')
    options = optimoptions(@quadprog,'Display','off');
    x = quadprog(HTH,f,[],[],[],[],lb,ub,[],options);
else
    func = @(x)0.5*norm(H*x-dw)^2;
    options.verbose = 0;
    options.method = 'newton';
    x = minConf_TMP(func,x0,lb,ub,options);
end

res=A*x - u;
wPrev = w;
wStep = wPrev + (r/m2)*(res - wPrev);

idxPos = wStep > 0;
if sum(idxPos) > k
    wPos = wStep(idxPos);
    [~,idxSort] = sort(wPos,'descend');
    wPos(idxSort(k+1:end)) = 0;
    wStep(idxPos) = wPos;
end
w=wStep;
if iter>maxIter
    break;
end
end


