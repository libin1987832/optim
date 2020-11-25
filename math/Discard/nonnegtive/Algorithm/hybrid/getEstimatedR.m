function rkn=getEstimatedR(QQ,r,I)
[m,n]=size(r);
rkn=r;
p=find(rkn>0);
 % r=(I-QQN)r
for i=1:m
    rkn(i)=r(i)-QQ(i,p)*r(p);
end
