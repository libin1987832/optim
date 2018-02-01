% 计算非负最小函数值 f=1/2||(b-Ax)+||
function f=fqf(b,A,x)
y0=b-A*x;
y0(y0<0)=0;
f=0.5*(y0'*y0);
