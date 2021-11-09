m = 10;
n = 5;

A = 2 * rand(m , n)-1;
colunmnormA=sum(A.^2,1);

pro1=colunmnormA/sum(colunmnormA)
residualvector=-b;
cumsumpro=cumsum(pro1)';
num_rand = 10000;
rr=rand(1,num_rand);
repmat(cumsumpro,1,n);
repmat(rr,n,1);
 l1=sum(repmat(cumsumpro,1,num_rand)<repmat(rr,n,1),1)+1;   %%%%% ¸ÅÂÊÑ¡È¡
m1=tabulate(l1);

index=1:n;
weight = colunmnormA/sum(colunmnormA);
picked = zeros(num_rand,1);
for i=1:num_rand
picked(i)=randsample(index,1,true,weight);
end;
picked';
m2=tabulate(picked);

iter_index=1;
for i=1:n-1
    k = ceil(weight(i)*num_rand);
    pickedj_a(iter_index:iter_index+k-1)=repmat(i,1,k);
    iter_index=iter_index+k;
end
pickedj_a(iter_index:num_rand)=repmat(n,1,num_rand-iter_index+1);
pickedj_i=randperm(num_rand);
m3=tabulate(pickedj_a(pickedj_i));
[pro1'*100,m1(:,3),m2(:,3),m3(:,3)]