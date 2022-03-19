% & CHA & 36.4281 & 1.29581e-12 & (11,1) & 11.831 \\
% & RHA & 36.4281 & 1.21063e-12 & (23,1) & 14.563 \\
% & PHA & 36.4281 & 1.29597e-12 & (9,1) & 11.924 \\
% & CCD & 36.4281 & 1.0207e-11 & 60000 & 2.59 \\
% & UCD & 36.4281 & 1.1723e-11 & 120000 & 5.299 \\
% & RCD & 36.4281 & 1.30059e-11 & 120000 & 5.637 \\

% & CHA & 27.8561 & 5.95211e-13 & (9,1) & 2.866 \\
% & RHA & 27.8561 & 5.74519e-13 & (23,1) & 4.164 \\
% & PHA & 27.8561 & 5.95211e-13 & (9,1) & 3.162 \\
% & CCD & 27.8561 & 5.25391e-12 & 40000 & 1.177 \\
% & UCD & 27.8561 & 9.94412e-12 & 100000 & 3.098 \\
% & RCD & 27.8561 & 1.00268e-11 & 100000 & 3.272 \\

% & CHA & 9.00499 & 1.69799e-14 & (36,1) & 0.015 \\
% & RHA & 9.00499 & 1.6933e-14 & (5,1) & 0.007 \\
% & PHA & 9.00499 & 1.81286e-14 & (2,1) & 0.008 \\
% & CCD & 9.00499 & 1.49994e-13 & 4000 & 0.013 \\
% & UCD & 9.00499 & 2.88041e-13 & 10000 & 0.05 \\
% & RCD & 9.00499 & 2.88945e-13 & 10000 & 0.039 \\
X = categorical({'Small 500X100','Medium 6000X600','Large 10^4X10^3'});
X = reordercats(X,{'Small 500X100','Medium 6000X600','Large 10^4X10^3'});
x=[500 6000 10000]
y1=([0.015,0.007,0.008,0.013,0.05,0.039])*30;
y2=([2.866,4.164,3.162,1.177,3.098,3.272]);
y3=([11.831,14.563,11.924 2.59,5.299,5.637]);

y=[y1;y2;y3];
% plot(x,y)
b=bar(X,y);
 ylabel('CPU运行时间');
 xlabel('矩阵的维数');
 Legend=legend('CHA','RHA','PHA','CCD','UCD','RCD');
title('在不同规模的矩阵下混合算法和坐标下降算法性能比较')
% for i =1 :6
% xtips2 = b(i).XEndPoints;
% ytips2 = b(i).YEndPoints;
% labels2 = string(b(i).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% end
% figure(10);%VAF Plot 
% VAFBar=rand(1,20);
% Bar=bar(VAFBar*100);
% BarH = findobj(Bar,'type','patch');
% hatchfill(BarH(1),'single',45,1,'r');
% hatchfill(BarH(2),'single',180,0.2,'b');
% hatchfill(BarH(3),'single',135,1,'g');  
% hatchfill(BarH(4),'cross',180,1,'m'); 
% hatchfill(BarH(5),'cross',135,1,'c'); 
% grid on;
% %ch = get(Bar,'children');
% set(gca,'XTickLabel',{'Uniform1','Vortex1','Uniform2','Vortex2'})
% % set(ch,'FaceVertexCData',[1 0 1;0 0 0;])
% %legend('Optimized A D','Optimized B D1','Optimized B D2','Optimized B D3','Optimized B D4','Optimized B D5','Location','SouthEast');
% Legend=legend('Optimized A D1','Optimized A D2','Optimized B D1','Optimized B D2','Optimized B D3',0);
% LegendH = findobj(Legend,'type','patch');
% hatchfill(LegendH(1),'single',45,1,'r');
% hatchfill(LegendH(2),'single',180,0.2,'b');
% hatchfill(LegendH(3),'single',135,1,'g');  
% hatchfill(LegendH(4),'cross',180,1,'m'); 
% hatchfill(LegendH(5),'cross',135,1,'c'); 
% xlabel('Conditions ','FontSize',26,'Fontname','times new Roman');
% ylabel('R^2 ( % )','FontSize',26,'Fontname','times new Roman');
% axis([0.5 4.5 95 100]);
% set(gca,'Fontname','times new Roman','FontSize',26);