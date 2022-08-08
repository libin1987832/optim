function out = FM_DeblurTV(f,Img,K,sigma,opts)

lam = opts.lam; 
rho = opts.rho; 
tol = opts.tol; 
Nit = opts.Nit;
nosnr = opts.Nosnr;
%theta= 0.09;

relError        = zeros(Nit,1);
psnrError      =relError;
ssimError      =relError;
[row, col]  = size(f);
u           = f;

eigK        = psf2otf(K,[row col]); %In the fourier domain
eigKtK      = abs(eigK).^2;
eigDtD      = abs(psf2otf(1 , [row, col])).^2;
KD=eigKtK+eigDtD; 
    for k = 1:Nit
        relKu=imfilter(u,K,'circular');
        relr1=f-sigma-relKu;
        relr1(relr1<0)=0;
        relr2=-f+relKu;
        relr2(relr2<0)=0;
        relr3=-u;
        relr3(relr3<0)=0;
        relr4=-255+u;
        relr4(relr4<0)=0;
        relError(k)= norm(relr1)^2+norm(relr2)^2+norm(relr3)^2+norm(relr4)^2;   
        if nosnr
        psnr_fun=psnr(u,double(Img));
        ssim_index=ssim(u,double(Img));
        psnrError(k)=psnr_fun;
        ssimError(k)=ssim_index(1);
        end
        u_old=u;
       % Ku=imfilter(u,K,'circular');
        Ku=relKu; 
        r1=f-sigma-Ku;
        r1(r1<0)=0;
        r2=-f+Ku;
        r2(r2<0)=0;
        r3=-u;
        r3(r3<0)=0;
        r4=-255+u;
        r4(r4<0)=0;
     
        r1f = imfilter(r1,-K,'circular');
         r2f = imfilter(r2,K,'circular');
        r3f = -r3;
        r4f = r4;
%         r3f = imfilter(r3,-1,'circular');
%         r4f = imfilter(r4,1,'circular');
        r14f=r1f+r2f+r3f+r4f;
        du=fft2(-r14f)./KD;
        du = real(ifft2(du));
        u=u_old+du;
%         relKu=imfilter(u,K,'circular');
%         relr1=f-sigma-relKu;
%         relr1(relr1<0)=0;
%         relr2=-f+relKu;
%         relr2(relr2<0)=0;
%         relr3=-u;
%         relr3(relr3<0)=0;
%         relr4=-255+u;
%         relr4(relr4<0)=0;
%         
      
         
    end
    
    out.sol                 = u;
    out.relativeError       = relError(1:k);
    out.psnrError = psnrError;
    out.ssimError = ssimError;
end


