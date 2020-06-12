clc
clear
addpath('./other')
addpath('./symsub')
addpath('../other')
addpath('../symsub')

n=2000;
C=sprandsym(n,0.4,1)+speye(n)*1.5;
xs=sprandn(n,1,0.3);
% xs must be nonnegative
I1=(xs<0);
xs(I1)=0;
% count the number of active set 
I2=(xs>0);
cs2=sum(I2);
% generate q
q=sprandn(n,1,0.3);
qt=C*xs;
q(xs>0)=-qt(xs>0);
% nondegenency
q(xs==0)=max(abs(qt))+0.1;

save('fpi','C','xs','q','n');
% load('fpiq')

% % convergence of the iteratve methods
% CC=triu(C,1);
% TC=tril(C);
% [symmetric2,posdef2]=isPosdef(TC-CC);
% if posdef2==0
% disp('the iterative methods may not converge');
% else
% disp('the iterative methods converges');
% end
% % object value
% ressvvf=func(C,q,xs);
% % KKT error
% ressvv=test_valid(C,q,xs);
% % posdef
% [symmetric,posdef]=isPosdef(C);
% if posdef==0
% disp('C may not the prositive semidefine');
% else if posdef==1
%     disp('positive');
% else
%     disp('semipositive');
%     end
% end
 x0=sparse(1:n,1:1,ones(n,1));
% x0=speye(n,1);
nmax=10;

max_iter = nmax;
tol_rel  = 1e-5;
tol_abs  = 1e-10;
nf=10;

[xkpsor err iter flag convergence msg] = psor(C, q, x0, 1, max_iter*nf/2, tol_rel, tol_abs, true);
[xks,ress]=splitS(C,q,1.4,x0,10);
disp(['mindig:' num2str(min(diag(C)))]);
[xk2,err,index2]=splitForlcp(x0,nmax,nf,C,q);
[xkpa,errpa,indexpa1,indexpa2]=PA(x0,nmax,nf,C,q);
[xkor,error,indexor,indexNor]=hybridorigin(x0,nmax,nf,C,q);
% || x-x^*||
xspsor=norm(xs-xkpsor);
xs2=norm(xs-xk2);
xspa=norm(xs-xkpa);
xsor=norm(xs-xkor);
% ||KKT||
[ress,fxs]=test_valid(C,q,xs);
[resss,fxss]=test_valid(C,q,xks);
[res2,fx2]=test_valid(C,q,xkpsor);
[res3,fx3]=test_valid(C,q,xk2);
[res4,fx4]=test_valid(C,q,xkpa);
[resor,fxor]=test_valid(C,q,xkor);

disp(['RFNP:(',num2str(index2*(nf+1)),' ',num2str(index2),') err:',num2str(res3)])
disp(['PA:(',num2str(indexpa1),' ',num2str(indexpa2),') err:',num2str(res4)])
disp(['our:(',num2str(indexor),' ',num2str(indexNor),') err:',num2str(resor)])
disp(['GS:(',num2str(max_iter),' ',num2str(0),') err:',num2str(res2)])


% % analysis
% checkEqS(cell2mat(error(4,1)),xs)
% checkEqS(cell2mat(error(4,2)),xs)
% checkEqS(cell2mat(error(4,3)),xs)
% checkEqS(cell2mat(error(4,4)),xs)
