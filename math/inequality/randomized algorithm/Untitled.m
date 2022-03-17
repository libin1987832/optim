%close all 
clc
clear

%parameters
m = 10;
n = 2; 
sigma = 0.1; 
p = 0.5; %Bernoulli probability  
numRuns = 10; %number of trials 
error_matrix = zeros(1000,numRuns); %preallocating error_matrix
error2_matrix = zeros(1000,numRuns); %preallocating error_matrix


for k = 1:numRuns
    %set up for matrix A 
    A = randn(m,n);
%     A(1000,:) = zeros(1,n); %uncomment to do inconsistent 
    b = randn(n,1); 
    y = A*b;  
   y = y + sigma*randn(m,1);
    b = A\y;
    
    %Bernoulli 
    B = (rand(m,n)<p); 
    
    %set up for matrix M
    M = A.*B; 
    
    %initialize x_0 and r_0
    x = zeros(n,1); %for Matrix A
    r = y - A*x; 
    x2 = zeros(n,1); %for Matrix M
    r2 = y - M*x;
    
    for t = 1:1000 %maximum number of iterations
        %Matrix A 
        j = datasample(1:n,1); %datasample(data,k) returns k observations sampled uniformly at random, with replacement, from the data 
        alpha_t =  (A(:,j)'*r)/ norm(A(:,j))^2; %' transposes vector/matrix
        e = zeros(n,1); 
        e(j) = 1; %jth euclidean basis vector 
        x = x + alpha_t * e; 
        r = r - alpha_t * A(:,j); 

        %Matrix M
        j2 = datasample(1:n,1);
        alpha_t2 = (M(:,j2)'*r2)/ norm(M(:,j2))^2;
        e2 = zeros(n,1);
        e2(j2) = 1; 
        x2 = x2 + alpha_t2 * e2; 
        r2 = r2-alpha_t2 * M(:,j2);
        
        error(t) = norm(b-x); %Matrix A
        error2(t) = norm(b - x2); %Matrix M
        
    end
    error_matrix(:,k) = error; 
    error2_matrix(:,k) = error2;
end 

%take a column-wise average to get the average error at each iteration
avgErr2 = mean(error2_matrix,2);
avgErr = mean(error_matrix,2);

%output
semilogy(avgErr, 'r') % Matrix A
xlabel('Iterations') 
ylabel('Error') 
hold on 
semilogy(avgErr2, 'b') % Matrix M
legend('Consistent Matrix A', ' Consistent Matrix M ', 'Inconsistent Matrix A','Inconsistent Matrix M' )
% MSE = (sum((x - x2).^2))/n 
% saveas(gcf, '5000 trials, 0.9p 0.9sigma.jpg')
% % plot (error , 'b')

%color code 
%Consistent Matrix A = red
%Consistent Matrix M = blue
%Inconsistent Matrix A = green
%Inconsistent Matrix M = black 
