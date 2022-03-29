
% Created by 12/12/2017 by Tarmizi Adam.
% This demo is to show ADMM total variation deblurring. The main solver is
% "ADMM_DeblurTV.m". Refer to mentioned file for more details

clc;
clear all;

ImgA={'Lena','male','mandril'};
name=['1','2','3'];
load('dd')
for i=1:3

set(gcf,'position',[100,100, 700, 600]);
set(gca,'DataAspectRatio',[1,1,1]); 
    % Img = imread([name(i) '.tiff']);
% 
% if size(Img,3) > 1
%     Img = rgb2gray(Img);
% end
% % ImgA{i}=Img;
% 
% K     =   fspecial('average',3); % For denoising
% f = imfilter(Img,K,'circular');
% 
% f = double(f);
% 
% BSNR = 20;
% sigma = BSNR2WGNsigma(f, BSNR);
% 
% f = f +  sigma * rand(size(Img)); %Add a little noise
% 
% %*** ADMM algorithm parameter set up ***
% 
% opts.lam = 1.5;
% opts.rho = 1.3;
% opts.tol = 1e-5;
% opts.Nit = 30;
% opts.Nosnr=0;
% out1 = FM_DeblurTV(f,Img,K,sigma,opts);
% out2 = ADMM_DeblurTV(f,Img,K,opts);
%subplot(4,3,i)
sigma=1;
x=(i-1)*0.2;
subplot('position',[0+x,0,.2,.2]);
imshow(Img)
title(sprintf('���������ŵ�%s',ImgA{i}))
% subplot(4,3,i+3)
subplot('position',[0+x,0.25,.2,.2]);
imshow(uint8(f))
title(sprintf('�������\\sigma=%g',sigma),'Interpreter','tex')
% % subplot(4,3,i+3)
subplot('position',[0+x,0.5,.2,.2]);
imshow(uint8(out2.sol))
title('ȫ��ֵ�ͼ���ع�')
% % subplot(4,3,i+3)
subplot('position',[0+x,0.75,.2,.2]);
imshow(uint8(out1.sol))
title('ȫ��ֵ�ͼ���ع�')
end                   % plot(out.relativeError)