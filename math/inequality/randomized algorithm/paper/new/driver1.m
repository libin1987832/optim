%driver 

%% experiment 1 (presentation) comparing rate of convergence among 3 versions with maxit = 15000 (same experiment as from the original paper)
maxit = 200;

% m = 150;
% n = 10;
r = 5;
m = 100;
n = 2*r+1;

% generate nodes
%tj are drawing randomly from a uniform distribution in [0,1]
t = rand(m,1);

%then ordering by magnitude (since they are all positive, just sort)
t = sort(t);
%just to assign the size
w = zeros(m,1); 
x = zeros(n,1);
% generate x 
realx = randn(n,1);
imgx = randn(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
% generate A and b
A = zeros(m,n);
for j = 1:m
%    dealing with special cases when reach the endpoints(nodes)
    if j == 1
        w(j) = (t(2)-t(end)-1)/2;
    elseif j == m
            w(j) = (t(1)+1-t(m-1))/2;
        else
            w(j) = (t(j+1)-t(j-1))/2;
    end              
    for k = -r:r
       A(j,k+r+1) =sqrt(w(j))*exp(2*pi*1i*k*t(j));
    end
end  
 A = 2 * rand(m , n)-1;
%A=real(A);
b = 2 * rand(m , 1)-1;

x0 = zeros(n , 1);

 x_exact=[];
nf=10;
maxIter = 2000;
tol=1e-20;
%[x_HA,flag,iter_HA,error_k,indexsm] = hybridA(A,b,x0,maxIter,nf,'PHA',tol,x_exact,0);
[x_HA, ~, ~, ~, ~] = IFM(A, b, x0,maxIter, nf , 1e-10,[],0);
r = b - A * x_HA;
r(r<0) = 0;
r_GS = r'*r;
g_GS = norm(A'*r)


figure (1)
iter=50;
maxit=390;
errors=zeros(3,iter);
% title('comparing rate of convergence among 3 versions with maxit = 15000')

for i = 1:iter
[error1,error2,error3] = driver_as_function(A,b,x0,r_GS,maxit,[]);

semilogy(1:maxit+1,error1/2,'k--','LineWidth',3);
% title('comparing rate of convergence among 3 versions with maxit = 15000')
% ylabel('Least squares error') 
% xlabel('Number of iterations') 
hold on
% semilogy(1:maxit+1,error2/2,'b--');
semilogy(1:maxit+1,error3/2,'r--');
% legend('classical kaczmarz','simple randomized kaczmarz','randomized kaczmarz')
errors(:,i)=[error1(1,maxit);error2(1,maxit);error3(1,maxit)];
end
 ylabel('F(x^k)-F^*') 
 xlabel('迭代次数') 
legend('循环模式','依列范数比例的概率分布')
hold off
% figure(2)
% hist(errors(1,:),10)
% figure(3)
% hist(errors(2,:),10)
% figure(4)
% hist(errors(3,:),10)
% %% experiment 2 (report) Running 100 times with random input and maxit = 5000 - statistical analysis
% error_c = [];
% error_s = [];
% error_r = [];
% for n = 1:100
%     [error1,error2,error3] = driver_as_function(5000,[]);
%     error_c = [error_c,error1(5000)];
%     error_s = [error_s,error2(5000)];
%     error_r = [error_r,error3(5000)];
% end
% % histogram
% figure (2)
% title('errors')
% h3 = histogram(error_c);
% ylabel('frequency') 
% xlabel('Least squares error at the 5000th iteration')
% hold on
% h4 = histogram(error_s);
% h5 = histogram(error_r);
% hold off
% legend('classical kaczmarz','simple randomized kaczmarz','randomized kaczmarz')
% 
% 
% figure (3)
% h1 = histogram(error_s);
% ylabel('Frequency') 
% xlabel('Least squares error at the 5000th iteration')
% hold on
% h2 = histogram(error_r);
% hold off
% legend('simple randomized kaczmarz','randomized kaczmarz')
% 
% %% experiment 3 (report) Running 100 times with random input and tol = 10^(-4) - statistical analysis
% iter_c = [];
% iter_s = [];
% iter_r = [];
% for n = 1:100
%     [iter1,iter2,iter3] = driver_as_function([],0.0001);
%     iter_c = [iter_c,iter1];
%     iter_s = [iter_s,iter2];
%     iter_r = [iter_r,iter3];
% end  
% % histogram
% figure (4)
% title('iterations')
% h6 = histogram(iter_c);
% ylabel('Frequency') 
% xlabel('Number of iterations to reach the tolerrance')
% hold on
% h7 = histogram(iter_s);
% h8 = histogram(iter_r);
% hold off
% legend('classical kaczmarz','simple randomized kaczmarz','randomized kaczmarz')
% figure (5)
% title('iterations')
% h9 = histogram(iter_s);
% ylabel('Frequency') 
% xlabel('Number of iterations to reach the tolerrance')
% hold on
% h10 = histogram(iter_r);
% hold off
% legend('simple randomized kaczmarz','randomized kaczmarz')