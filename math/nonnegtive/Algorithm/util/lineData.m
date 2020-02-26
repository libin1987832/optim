function d=lineData(A,b)
[m,n]=size(A);
d=zeros(m,4);
for i=1:m
    if abs(A(i,1))>1e-8
        d(i,1)=b(i,1)/A(i,1);
    else
        d(i,1)=1;
        if abs(A(i,2))>1e-8
            d(i,2)=b(i,1)/A(i,2);
        else
            print('zores')
        end
    end
    if abs(A(i,2))>1e-8
        d(i,4)=b(i,1)/A(i,2);
    else
        d(i,3)=1;
        if abs(A(i,2))>1e-8
            d(i,4)=b(i,1)/A(i,1);
        else
            print('zores')
        end
    end
end