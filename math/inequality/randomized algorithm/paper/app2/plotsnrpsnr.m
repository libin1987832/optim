
% Created by 12/12/2017 by Tarmizi Adam.
% This demo is to show ADMM total variation deblurring. The main solver is
% "ADMM_DeblurTV.m". Refer to mentioned file for more details

clc;
clear all;
close all;

Img = imread('1.tiff');

if size(Img,3) > 1
    Img = rgb2gray(Img);
end

K     =   fspecial('average',5); % For denoising
f = imfilter(Img,K,'circular');

f = double(f);

BSNR = 20;
sigma = BSNR2WGNsigma(f, BSNR)*2;
f = f +  sigma * rand(size(Img));
imshow(uint8(f));
% u=Img;
% relKu=imfilter(u,K,'circular');
% relKu=double(relKu);
% relr1=f-sigma-relKu;
% relr1(relr1<0)=0;
% relr2=-f+relKu;
% relr2(relr2<0)=0;
% relr3=double(-u);
% relr3(relr3<0)=0;
% relr4=double(-255+u);
% relr4(relr4<0)=0;
% relError= norm(relr1)^2+norm(relr2)^2+norm(relr3)^2+norm(relr4)^2
% f = f +  sigma * randn(size(Img)); %Add a little noise


%*** ADMM algorithm parameter set up ***

opts.lam = 1.5;
opts.rho = 1.3;
opts.tol = 1e-5;
opts.Nit = 20;
opts.Nosnr = 1;
out1 = FM_DeblurTV(f,Img,K,sigma,opts);
out2 = ADMM_DeblurTV(f,Img,K,opts);
 figure;
% semilogy(out1.relativeError)
% figure
%  semilogy(out2.relativeError)
% out1.relativeError(end)
plot(1:opts.Nit,out1.psnrError,1:opts.Nit,out2.psnrError );
xlabel('迭代次数')
ylabel('峰值信噪比(PSNR)')
legend('最小二乘意义下线性不等式方程组','全变分模型')
title('最小二乘意义下线性不等式方程组重构Lena图像');
figure;
title('最小二乘意义下线性不等式方程组重构Lena图像');
plot(1:opts.Nit,out1.ssimError,1:opts.Nit,out2.ssimError);
xlabel('迭代次数')
ylabel('结构相似性(SSIM)')
legend('最小二乘意义下线性不等式方程组','全变分模型')
% figure
% imshow(uint8(out1.sol))
% figure
% imshow(uint8(out2.sol))
% psnr_fun=psnr(out1.sol,double(Img));
% ssim_index=ssim(out1.sol,double(Img));

% title(sprintf('ADMM-TV Deblurred (PSNR = %3.3f dB,SSIM = %3.3f) ',...
%                        psnr_fun,ssim_index(1)));
%                    figure()
 %                  plot(out1.relativeError)

                   % plot(out.relativeError)