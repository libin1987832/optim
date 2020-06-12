function [alph,xk,res]=projectedsearch(d,x0,G,c)
tol=1e-10;
res=0;
xk=x0;
alph=0;
if norm(d)>tol
    [m,n]=size(G);
    le0=(d<0);
    % negative is marked for inf 
    t1=-1*ones(n,1);
    % breakpoints
    t1(le0)=x0(le0)./(-d(le0));
    % d>0 d==0
    maxt1=max(t1);
    % inf
    t1(t1<0)=maxt1+1;
    st=sort(unique([t1;0]));
    len=length(st)-1;
    %[i-1,i] trie
    for i=1:len
        t=t1;
        % d all step length
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
    res=1;
    % test_projected(alph,G,c,x0,d)
end
