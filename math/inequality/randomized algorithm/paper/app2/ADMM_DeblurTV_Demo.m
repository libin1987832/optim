
% Created by 12/12/2017 by Tarmizi Adam.
% This demo is to show ADMM total variation deblurring. The main solver is
% "ADMM_DeblurTV.m". Refer to mentioned file for more details

clc;
clear all;

ImgA={'Lena','male','mandril'};
name=['1','2','3'];

for i=1:3

set(gcf,'position',[100,100, 700, 600]);

    Img = imread([name(i) '.tiff']);

if size(Img,3) > 1
    Img = rgb2gray(Img);
end
% ImgA{i}=Img;

K     =   fspecial('average',20); % For denoising
f = imfilter(Img,K,'circular');

f = double(f);

BSNR = 20;
sigma = BSNR2WGNsigma(f, BSNR);

f = f +  sigma * rand(size(Img)); %Add a little noise

%*** ADMM algorithm parameter set up ***

opts.lam = 1.5;
opts.rho = 1.3;
opts.tol = 1e-5;
opts.Nit = 100;
opts.Nosnr=0;
out1 = FM_DeblurTV(f,Img,K,sigma,opts);
out2 = ADMM_DeblurTV(f,Img,K,opts);
%subplot(4,3,i)

x=(i-1)*0.2;
subplot('position',[0+x,0.75,.2,.2]);
% subfig=subplot(4,3,i);
% RemoveFigMarginArea(subfig)
% RemoveSubplotMarginArea(subfig, 4,3,i)
imshow(Img)
% set(gca,'DataAspectRatio',[3,4 ,1]);  
title(sprintf('无噪声干扰的%s',ImgA{i}))
%subfig=subplot(4,3,i+3)
%RemoveFigMarginArea(subfig)
% RemoveSubplotMarginArea(subfig, 4,3,i+3)
subplot('position',[0+x,0.5,.2,.2]);
% set(gca,'DataAspectRatio',[3,4,1]);  
imshow(uint8(f))
title(sprintf('添加噪声\\sigma=%g',sigma),'Interpreter','tex')
% subfig=subplot(4,3,i+6)
% RemoveFigMarginArea(subfig)
%  RemoveSubplotMarginArea(subfig, 4,3,i+6)
subplot('position',[0+x,0.25,.2,.2]);
% set(gca,'DataAspectRatio',[3,4,1]);  
imshow(uint8(out2.sol))
title('全变分图像重构') 
% subfig=subplot(4,3,i+9)
%  RemoveFigMarginArea(subfig)
% set(subfig,'LooseInset',get(subfig,'TightInset'))
% RemoveSubplotMarginArea(subfig, 4,3,i+9)
subplot('position',[0+x,0,.2,.2]);
% set(gca,'DataAspectRatio',[3,4,1]);  
imshow(uint8(out1.sol))
title('不等式图像重构')

end                   % plot(out.relativeError)

fdata = getframe(gcf);
figure
imshow(fdata.cdata(:,2:425,:))
%imwrite(f.cdata(:,10:420,:), 'imagediff.fig');
%imwrite(f.cdata(:,10:420,:), 'imagediff.eps');
