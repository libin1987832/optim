%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
addpath('./dataInequality/');
addpath('./algorithmInequality/');
addpath('./randomized algorithm')

m = 1000;
n = 20;
rangeMax = 2;
rangeMin = -2;

   A = 2 * rand(m , n)-1;
   b = 2 * rand(m , 1)-1;
  % b=A*ones(n,1);
   x0 = zeros(n , 1);
   maxIter = 300;
   nf = 5;
   str = ['D','U','C','R','P'];
% for the solution so here

    [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
    [rk, rkh, dh, gh] = residual(A,b,xkh);
    fprintf('active:%d ,%g',vkh,dh);
    maxit = 10000;
    [xkacz,iterkacz,errorkacz] = randomizedKaczmarzNE(A, b, x0,2*maxit,[],xkh);
    [xGS,iterGS,errorGS] = randomizedGaussSeidelNE(A, b, x0,maxit,[],xkh);



 %%

beginp = 1;
figure
%maxIterA = 200;
h=semilogy(beginp:iterkacz,errorkacz,'b.');
h.LineStyle = '--';
hold on
h=semilogy(beginp:iterGS,errorGS);
h.LineStyle = '--';
legend('Kaczmarz','Gauss Seidel');
xlabel('the iterative numbers');
ylabel('the norm of the error');

%    end % for ration
%    end % for m