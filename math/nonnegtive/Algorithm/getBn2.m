function sum2=getBn2(nIter,Q,r,I)
[m,n]=size(r);
rkn=r;
%p=find(rkn>0);
% sum=0;
sum2=0;
for i=1:m
    tmp=0;
    for j=1:m
        if r(j)>0
            tmp=tmp+Q(i,:)*Q(j,:)'*r(j);
        end
    end
    rkn(i)=r(i)-nIter*tmp;
%     sum=sum+sign(rkn(i)*r(i));
    if sign(rkn(i)*r(i))>0
        sum2=sum2+1;
    end
%     if sum2~=sum
%         error('dddd')
%     end
end
