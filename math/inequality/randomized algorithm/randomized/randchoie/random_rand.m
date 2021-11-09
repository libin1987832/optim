function pron = random_rand(pro,num_rand)
n=size(pro,2);
cumsumpro=cumsum(pro);
pron=sum(repmat(cumsumpro,num_rand,1)<repmat(rand(num_rand,1),1,n),2)+1;   %%%%% ¸ÅÂÊÑ¡È¡