function sum2=getABn(QQ,r,I,niter)
% sum2 is to conuting number which is same r
% QQ=Q'*Q
% Nk=r;
% Nk(Nk>0)=1;
% Nk(Nk<0)=0;
% Bn=I-QQ*diag(Nk);
% rkn2=Bn*r;
[m,n]=size(r);
rkn=r;
IQQN=I-QQ/niter*diag(rkn>0);
IQQNn=IQQN^(50*niter);
sum=0;
sum2=0;
 % r=(I-QQN)r
for i=1:m
    rkn(i)=IQQNn(i,:)*r;
%     sum=sum+sign(rkn(i)*r(i));
    if sign(rkn(i)*r(i))>0
        sum2=sum2+1;
    end
end
