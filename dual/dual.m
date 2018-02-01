% 对偶算法 算法见：A numerically stable dual method for solving strictly convex quadratic problems
invG=inv(G)
x=-invG*a;
f=0.5*a'*x;
H=invG;
A=[];
u=[];
q=0;
S=C'*x-b;
Nx=[];
[V,n]=find(S<0);
if sum(V)>0
p=V(1);
nj=C(p);
uj=[u;0];
if 0==q
   u=0
end
z=H*nj;
if q>0 
    r=Nx*nj;
end
if q==0 || r<=0
    t1=+inf;
else
    t1=+inf;
    col=0;
    for i=1:n
        if r(i)>0
            tt=uj(i)/r(i);
            if tt<t1
                t1=tt;
                col=i;
            end
        end
    end
end
if norm(z)==0 
    t2=+inf;
else
    t2=-S(p)/z'*nj;
end
t=min(t1,t2);

        
        
