% n=10000;
% tic
% A=sprandsym(n,0.01,1e-10,2);
% f=toc;
% % condest(A);
% fprintf('%8.04f\n',f)
% xs=sprandn(n,1,0.3);
% b=A(5,:)*xs;
% fa=full(A(5,:))*full(xs);
% tic
% b=A*xs;
% f1=toc;
% Ic=cell(n,1);
% for i=1:n
%     Ic(i,1)={find(A(i,:)>-1e-10)};
% end
% c=zeros(n,1);
% tic
% for i=1:n
%     c(i)=A(i,cell2mat(Ic(i,1)))*xs(cell2mat(Ic(i,1)));
% end
%     
% f2=toc;
% sum(c-b)
% fprintf('%8.04f,%8.04f,%8.04f\n',f,f1,f2)
A = sprand(10000,10000,0.005);
B = sprand(10000,1,0.005);
% Af = full(A);
% Bf = full(B);
% timeit(@() Af*Bf)

timeit(@() A*B)
timeit(@() A(1,1:10000)*B)
