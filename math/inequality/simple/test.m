



 semilogx(tolsa,arrspeed(1,[10:-1:1]),'b*-',tolsa,arrspeed(2,[10:-1:1]),'ro-',tolsa,arrspeed(3,[10:-1:1]),'g+-') 
% semilogx(tolsa,arrspeed(1,:),'b*-',tolsa,arrspeed(2,:),'ro-',tolsa,arrspeed(3,:),'g+-')
  gca=axes('XDir','reverse');
    gca=axes();
 set(gca,'XTick',tolsa)
%  set(gca,'XDir','reverse')
% % 标题标注 
% title('Compared with other algorithms, the speed up of CHA is relation with the accuracy') 
% % 坐标轴标注 
% xlabel('the accuracy') 
% ylabel('the speed up') 
% legend('FM/CHA','IFM/CHA','HAN/CHA') 