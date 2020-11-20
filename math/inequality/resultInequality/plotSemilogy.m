function plotSemilogy(beginp,maxIterA,iterA,type,typet)
%% plot picture
type=['r','c','k','g'];
typet=['+','o','v','s'];
beginp = 1;
figure
%maxIterA = 70;
h=semilogy(beginp:maxIterA,iterA(1,beginp:maxIterA),'bx');
h.LineStyle = '--';
hold on
for i=2:4
    h=semilogy(beginp:maxIterA,iterA(i,beginp:maxIterA),[type(i) typet(i)]);
    h.LineStyle = '--';
end
legend('DHA','GHA','RHA','PHA');
xlabel('Iteration Number');
ylabel('the norm of the gradient');