

function SKMtest(A,b,numSolves,maxIter,error_bound,beta,lambda,sysTitle)


global data    %saves relevant information per test: 
               %beta, lambda, average time, average number of iterations
            
global ITER;   %saves the number of iterations to terminate
	           %the algorithm (reach relative residual error or MaxIter)		   

[m,n]=size(A);
k=size(beta,2); %number of sample sizes to be tested 
lamb = size(lambda,2); %number of projection constants to be tested

data = zeros(4,k,lamb); 

for i = 1:lamb               %for all the values of lambda
        for j=1:k         %for all the values of beta
             disp(['Running for beta = ',num2str(beta(j)),...
             'lambda = ',num2str(lambda(i))])
	         timer=[];
             ITER=[];
             for count = 1:numSolves
		         x_0 = 10*ones(n,1); 
              	 tstart=cputime;
                 SKM(A,b,x_0,lambda(i),error_bound,maxIter,beta(j));
                 timer = [timer,cputime-tstart];
             end
             
             %saves beta 
             data(1,j,i) = beta(j); 	       
             
             %saves lambda
             data(2,j,i) = lambda(i); 
             
             %saves mean of the number of iterations rounded down
             data(3,j,i) = floor(mean(ITER));    
             
             %saves mean of the computational time
	         data(4,j,i) = mean(timer); 	         

        end
       
end

    save([sysTitle,'Experimentwithm=', num2str(m), 'n=', num2str(n), '.mat']);

    %plot of computational time as a function of the sample size
    figure(1);
    
    semilogy(data(1,:,1), data(4,:,1), 'k', data(1,:,2), data(4,:,2), ...
    '-.b', data(1,:,3), data(4,:,3), 'r', data(1,:,4),data(4,:,4),'--g', ...
    data   (1,:,5),data(4,:,5),'m','LineWidth',3);
    set(gca,'fontsize',16);

    axis([1,inf,-inf,inf]);
    
    title(['Netlib LP ', sysTitle, ' feasibility problem: ', num2str(m),' x ',num2str(n)],'fontsize',20);

    xlabel('Sample Size, \beta','fontsize',18);
    ylabel(['Comp. Time for relative residual halting error: ', ...
    num2str(error_bound)],'fontsize',18);

    legend(['\lambda = ',num2str(lambda(1))], ...
    ['\lambda = ', num2str(lambda(2))], ['\lambda = ',num2str(lambda(3))],...
    ['\lambda = ',num2str(lambda(4))], ['\lambda = ',num2str(lambda(5))]);

    print(['-f',num2str(1)],[sysTitle,'computational_time'],'-dpng')
    
    %plot of number of iterations as a function of the sample size
    figure(2);
    semilogy(data(1,:,1), data(3,:,1), 'k', data(1,:,2), data(3,:,2), ...
    '-.b', data(1,:,3), data(3,:,3), 'r',data(1,:,4),data(3,:,4),'--g', ...
    data(1,:,5),data(3,:,5),'m','LineWidth',3);
    set(gca,'fontsize',16);

    axis([1,inf,-inf,inf]);
    
    title(['Netlib LP ', sysTitle, ' feasibility problem: ', num2str(m),' x ',num2str(n)],'fontsize',20);

    xlabel('Sample Size,\beta','fontsize',18);
    ylabel(['Iterations for relative residual halting error: ', ...
    num2str(error_bound)],'fontsize',20);

    legend(['\lambda = ',num2str(lambda(1))], ...
    ['\lambda = ',num2str(lambda(2))], ['\lambda = ',num2str(lambda(3))],...
    ['\lambda = ',num2str(lambda(4))], ['\lambda = ',num2str(lambda(5))]);

    print(['-f',num2str(2)],[sysTitle,'iterations'],'-dpng')

end
