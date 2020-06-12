
m=100;
corA1=zeros(3,6);
corA2=zeros(3,6);
for n=100:100:600
    for nf=6:8
        cor=test_predict(n,m,nf);
        corA1(nf-5:(n-100)/100+1)=cor(1)/m;
        corA2(nf-5:(n-100)/100+1)=cor(2)/m;
    end
end
plot(100:100:600,corA1(:,1),100:100:600,corA1(:,2),100:100:600,corA1(:,3))