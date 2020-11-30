M = [1,0,0;2,1,0;2,2,1];
q = [-8;-12;-14];
[w,z,retcode] = Lemka(M,q)

% M = [10,0,-2;2,0.1,-0.4;0,0.2,0.1];
% q = [10;1;-1];

[w2,z,retcode] = Murty(M,q);
[w3,z,retcode] = Graves(M,q);
[w1,z,retcode] = Bard(M,q);
[w1';w2';w3']