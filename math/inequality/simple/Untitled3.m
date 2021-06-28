figure
% plot([0 nfA],[tfH records(1,:,5)],'k.-',[0 nfA],[tfH records(2,:,5)],'b*-',[0 nfA],[tfH records(3,:,5)],'ro-',[0 nfA],[tfH records(4,:,5)],'g+-')
 plot(nfA,records(1,:,5),'k.-',nfA,records(2,:,5),'b*-',nfA,records(3,:,5),'ro-',nfA,records(4,:,5),'g+-')
% % set(gca,'XTickLabel',tolsa);
% % 标题标注
 set(gca,'YLim',[0.2  0.6 ]);%X轴的数据显示范围
set(gca,'ytick',[0.2 0.3 0.4 0.5 0.6 1.3]);
 title('The performance of hybrid algorithms with increasing n_f') 
% % 坐标轴标注 
xlabel('n_f') 
ylabel('CPU(s)') 
legend('DHA(μ= n_f)','CHA','RHA','PHA') 
