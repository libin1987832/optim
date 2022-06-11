%% ͼ��AΪ�Գƾ���B��λ����Rou����1
x=importdata('XZB1.txt');
y=importdata('YZB1.txt');
x1=linspace(min(x),max(x));
y1=interp1(x,y,x1,'pchip');
semilogy(x1,y1,'r-','LineWidth',1.5)%������zhi��dao��
hold on

x=importdata('XZB5.txt');
y=importdata('YZB5.txt');
x5=linspace(min(x),max(x));
y5=interp1(x,y,x5,'pchip');
semilogy(x5,y5,'c-','LineWidth',1.5)%������zhi��dao��
hold on

x=importdata('XZB2.txt');
y=importdata('YZB2.txt');
x2=linspace(min(x),max(x));
y2=interp1(x,y,x2,'pchip');
semilogy(x2,y2,'m-','LineWidth',1.5)%������zhi��dao��

hold on
x=importdata('XZB4.txt');
y=importdata('YZB4.txt');
x4=linspace(min(x),max(x));
y4=interp1(x,y,x4,'pchip');
semilogy(x4,y4,'g-','LineWidth',1.5)%������zhi��dao��

hold on 
x=importdata('XZB3.txt');
y=importdata('YZB3.txt');
x3=linspace(min(x),max(x));
y3=interp1(x,y,x3,'pchip');
semilogy(x3,y3,'b-','LineWidth',1.5)%������zhi��dao��

title('B=C^TC+I')
axis([0,500,10E-7,10E+1])
legend('BAS','SPL','SQP(G)','SSQP(D)','SSQP(I)')
xlabel('Iterations')
ylabel('norm(dk)')