function [x,B,A]=nt_yulewalk_whiten(x,order,freqs)
%[y,B,A]=nt_yulewalk_whiten(x,order,freqs) - whiten spectrally
%
%  y: whitened data
%  B,A: filter coeffs
%
%  x: data to whiten
%  order: order of Yule-Walker filter
%  freqs: target frequencies (1 = Nyquist)
%
%  If freqs is not specified, distribute in intervals of equal power

NFFT=1024;

% estimate power spectrum
[s,f]=nt_spect_plot(x,NFFT,[],[],1);  

f=linspace(0,1,numel(f))'; % make frequency array start at 0

% whitening filter
[B,A]=yulewalk(order,f,1./sqrt(s));
x=filter(B,A,x);
 