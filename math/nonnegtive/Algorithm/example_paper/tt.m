function y=tt(n)
y=0;
s=-1;
for i=1:n;
  s=-1*s;
  y=y+1/i*s;
end