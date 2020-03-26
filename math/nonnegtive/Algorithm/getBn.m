function sum=getBn(QQ,r,I)
% Nk=r;
% Nk(Nk>0)=1;
% Nk(Nk<0)=0;
% Bn=I-QQ*diag(Nk);
% rkn=Bn*r;
[m,n]=size(r);
rkn=r;
p=find(rkn>0);
sum=0;
for i=1:m
    rkn(i)=r(i)-QQ(i,p)*r(p);
    sum=sum+sign(rkn(i)*r(i));
%     if sign(rkn(i)*r(i))>0
%         sum=sum+1;
%     end
end
