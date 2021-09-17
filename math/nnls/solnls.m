function out = solnls( A, b,  c, opt )
funObj2 = @(x) fungun(A,b,x);
funProj2 = @(x) proj(x);
x =c;
if strcmp(opt.algo, 'PLB')
[x,f,funEvals] = minConF_PQN(funObj2,x,funProj2,opt);
end
if strcmp(opt.algo , 'BB')
[x,f,funEvals,projects] = minConF_SPG(funObj2,x,funProj2,opt);
end
out.x = x;