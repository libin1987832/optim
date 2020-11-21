function plotSemilogy(res,type,typet)
%% plot picture
[m,n]=size(res);
figure
%maxIterA = 70;
h=semilogy(1:n,res(1,1:n),'bx');
h.LineStyle = '--';
hold on
for i=2:m
    h=semilogy(1:n,res(i,1:n),[type(i) typet(i)]);
    h.LineStyle = '--';
end
xlabel('Iteration Number');
ylabel('the norm of the gradient');