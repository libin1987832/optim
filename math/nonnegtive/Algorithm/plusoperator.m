function a=plusoperator(A,b,h,x)
A=-1*A;
a=1;
[m,n]=size(A);
xc=A*h;
bc=b-A*x;
ar=bc./xc;
compareAr=unique(sort([ar;0;1]));
M=zeros(m,1);
for i=1:n+2
    if compareAr(i)<1.0000001 && compareAr(i)>0.0000001
        tar=logical(ar<compareAr(i));
        M(logical(xc(tar)>-0.0001))=1;
        tar=logical(ar>compareAr(i));
        M(logical(xc(tar)<0))=1;
        M=diag(M);
        at=(xc'*M*bc)/(xc'*M*xc);
        if at<compareAr(i) && at>compareAr(i-1)
            a=at;
            break;
        end
    end
end


%test eg A=[-1;1] b=[0;1] x=0 h=1 a=1/2 plusoperator([1;-1],[0;1],1,0)