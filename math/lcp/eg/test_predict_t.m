 
m=500;
ns=100;niv=2;ne=120;
nns=4;nne=12;
corA1=zeros(nne-nns,(ne-ns)/niv);
corA2=zeros(nne-nns,(ne-ns)/niv);
for n=ns:niv:ne 
    for nf=nns:nne
        cor=test_predict(n,m,nf);
        corA1(nf-nns+1,(n-ns)/niv+1)=cor(1)/m;
        corA2(nf-nns+1,(n-ns)/niv+1)=cor(2)/m;
    end
end
corA1
plot([ns:niv:ne],corA1(1,:),[ns:niv:ne],corA1(2,:),[ns:niv:ne],corA1(3,:))