% M = [1,0,0;2,1,0;2,2,1];
% q = [-8;-12;-14];
% [w,z,retcode] = LCPSolveL(M,q)

M = [10,0,-2;2,0.1,-0.4;0,0.2,0.1];
q = [10;1;-1];
% [w,z,retcode] = Bard(M,q)
% [w,z,retcode] = Murty(M,q)
[w,z,retcode] = Graves(M,q)