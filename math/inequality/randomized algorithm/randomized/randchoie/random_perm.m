function pron = random_perm(pro,num_rand)
iter_index=1;
for i=1:n-1
    k = ceil(pro(i)*num_rand);
    pickedj_a(iter_index:iter_index+k-1)=repmat(i,1,k);
    iter_index=iter_index+k;
end
pickedj_a(iter_index:num_rand)=repmat(n,1,num_rand-iter_index+1);
pickedj_i=randperm(num_rand);
pron = pickedj_a(pickedj_i);