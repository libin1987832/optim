%% 图像，A为对称矩阵，B单位矩阵，Rou等于1
x=importdata('XFB1.txt');
y=importdata('YFB1.txt');
x1=linspace(min(x),max(x));
y1=interp1(x,y,x1,'pchip');
semilogy(x1,y1,'r-','LineWidth',1.5)%处理后的zhi曲dao线

hold on 
x=importdata('XFB5.txt');
y=importdata('YFB5.txt');
x5=linspace(min(x),max(x));
y5=interp1(x,y,x5,'pchip');
semilogy(x5,y5,'c-','LineWidth',1.5)%处理后的zhi曲dao线

hold on
x=importdata('XFB2.txt');
y=importdata('YFB2.txt');
x2=linspace(min(x),max(x));
y2=interp1(x,y,x2,'pchip');
semilogy(x2,y2,'m-','LineWidth',1.5)%处理后的zhi曲dao线

hold on
x=importdata('XFB4.txt');
y=importdata('YFB4.txt');
x4=linspace(min(x),max(x));
y4=interp1(x,y,x4,'pchip');
semilogy(x4,y4,'g-','LineWidth',1.5)%处理后的zhi曲dao线

hold on 
x=importdata('XFB3.txt');
y=importdata('YFB3.txt');
x3=linspace(min(x),max(x));
y3=interp1(x,y,x3,'pchip');
semilogy(x3,y3,'b-','LineWidth',1.5)%处理后的zhi曲dao线

title('B=C^TC+I')
axis([0,500,10E-7,10E+1])
legend('BAS','SPL','SQP(G)','SSQP(D)','SSQP(I)')
xlabel('Iterations')
ylabel('norm(dk)')