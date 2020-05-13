function [alph,xk]=projectedsearch(d,x0,G,c)
[m,n]=size(G);
le0=(d<0);
t1=-1*ones(n,1);
t1(le0)=x0(le0)./(-d(le0));
maxt1=max(t1);
t1(t1<0)=maxt1+1;
st=sort(unique([t1;0]));
len=length(st)-1;
for i=1:len
    t=t1;
    t(t1>st(i))=st(i);
    xt=x0+t.*d;
    p=zeros(n,1);
    p(t1>st(i))=d(t1>st(i));
   % fj=c'*xt+0.5*xt'*G*xt;
    fj1=c'*p+xt'*G*p;
    fj2=p'*G*p;
    if fj1>0
        alph=st(i);
        break;
    else
        deltt=-fj1/fj2;
        if st(i)+deltt<st(i+1) || st(i+1)>maxt1
            alph=st(i)+deltt;
            break;
        end
    end
end
xk=max(x0+alph*d,zeros(n,1));
% test_projected(alph,G,c,x0,d)