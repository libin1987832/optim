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
p=find(rkn>0);
sum=0;
sum2=0;
 % r=(I-QQN)r
for i=1:m
    rkn(i)=r(i)-QQ(i,p)*r(p);
%     sum=sum+sign(rkn(i)*r(i));
    if sign(rkn(i)*r(i))>0
        sum2=sum2+1;
    end
%     if sum2~=sum
%         error('dddd')
%     end
end
