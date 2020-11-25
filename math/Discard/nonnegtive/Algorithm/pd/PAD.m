function [xk,fk]=PAD(x0,A,b,NE,k)
Aa=A(NE,:);
Ba=b(NE,:);
[uk,fk]=krylov(Aa,Ba,k);
xk=uk;
fk=b-A*xk;
fk(fk<0)=0;
fk=0.5*fk'*fk;