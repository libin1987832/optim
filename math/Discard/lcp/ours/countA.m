function n=countA(x)
  I=(x>0);
  nt=zeros(size(x));
  nt(x>0)=1;
  n=sum(nt);