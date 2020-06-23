function varargout=getBn2(nIter,Q,r,varargin)
[m,n]=size(r);
rkn=r;
%p=find(rkn>0);
% sum=0;
sum2=0;
vn=size(varargin,1);
if vn > 0
    tmpq=varargin{1};
end
for i=1:m
    tmp=0;
    for j=1:m
        if r(j)>0
            if vn>0 && tmpq(i,j)<2
                tmp=tmp+tmpq(i,j)*r(j);
            else
                tmpq(i,j)=Q(i,:)*Q(j,:)';
                tmp=tmp+tmpq(i,j)*r(j);
            end
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
if nargout == 1
    varargout{1}=sum2;
else
    varargout{1}=sum2;
    varargout{2}=tmpq;
end
