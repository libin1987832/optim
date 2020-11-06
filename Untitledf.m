for i =1 :1000
A=rand(4,4)*10+123;
ATA = A'*A;
l1 = max(eig(ATA));
B = rand(1,4);
B = round(B);
l2 = max(eig(ATA*diag(B)));
if l1 > l2
    save('ss','A','B')
    [l1 l2];
end
end

A =  [1 -3;-3 1];
C=A'*A*[1 0;0 0];
eig((C+C')/2)
