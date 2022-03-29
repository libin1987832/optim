
% Created by 12/12/2017 by Tarmizi Adam.
% This demo is to show ADMM total variation deblurring. The main solver is
% "ADMM_DeblurTV.m". Refer to mentioned file for more details

clc;
clear all;
close all;

Img = imread('Lena_std.tif');

if size(Img,3) > 1
    Img = rgb2gray(Img);
end

K     =   fspecial('average',3); % For denoising
f = imfilter(Img,K,'circular');

f = double(f);

BSNR = 20;
sigma = BSNR2WGNsigma(f, BSNR)

f = f +  sigma * randn(size(Img)); %Add a little noise

%*** ADMM algorithm parameter set up ***

opts.lam = 1.5;
opts.rho = 1.3;
opts.tol = 1e-5;
opts.Nit = 100;

out1 = FM_DeblurTV(f,K,sigma/2,opts);
%out = ADMM_DeblurTV(f,Img,K,opts);
figure;
imshow(uint8(f));
psnr_fun=psnr(out.sol,double(Img));
ssim_index=ssim(out.sol,double(Img));
figure;
imshow(uint8(out.sol))
title(sprintf('ADMM-TV Deblurred (PSNR = %3.3f dB,SSIM = %3.3f) ',...
                       psnr_fun,ssim_index(1)));
figure;
psnr_fun=psnr(out1.sol,double(Img));
ssim_index=ssim(out1.sol,double(Img));
imshow(uint8(out1.sol))
title(sprintf('ADMM-TV Deblurred (PSNR = %3.3f dB,SSIM = %3.3f) ',...
                       psnr_fun,ssim_index(1)));

                   % plot(out.relativeError)