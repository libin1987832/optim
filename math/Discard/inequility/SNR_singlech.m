function snr=SNR_singlech(l,ln)
snr=0;
Ps=sum(sum((1-mean(mean(l))).^2));
Pn=sum(sum((1-ln).^2));
snr=10*log10(Ps/Pn);