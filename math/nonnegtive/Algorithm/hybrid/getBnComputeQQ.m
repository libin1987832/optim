% sum2 is to conuting number which is same r
% QQ=n*Q'*Q
function sum2=getBn(QQ,r,I)
% Nk=r;
% Nk(Nk>0)=1;
% Nk(Nk<0)=0;
% Bn=I-QQ*diag(Nk);
% rkn2=Bn*r;
[m,n]=size(r);
rkn=r;
% p can reduce computation only multiply nonzero
p=find(rkn>0);
sum=0;
sum2=0;
 % r=(I-QQN)r
for i=1:m
    rkn(i)=r(i)-QQ(i,p)*r(p);
    if sign(rkn(i)*r(i))>0
        sum2=sum2+1;
    end
end

