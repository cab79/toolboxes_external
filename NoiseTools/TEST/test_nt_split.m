clear
N=10000;

% single channel, step
x=[randn(N,1)+1; randn(N,1)-1]; 

figure(1); clf;
subplot 411
plot(x);
subplot 423
[idx,score_vector,score]=nt_split(x); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('find step')
subplot 424
[idx,score_vector,score]=nt_split(x.^2); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('find step in power')

% single channel, step in power
x=[randn(N,1); 2*randn(N,1)];  
subplot 413
plot(x);
subplot 427
[idx,score_vector,score]=nt_split(x); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('find step')
subplot 428
[idx,score_vector,score]=nt_split(x.^2); 
plot(score_vector); title('find step in power')
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   


% 10 channels
a=nt_normcol(randn(N,5)*randn(5,10));
b=nt_normcol(randn(N,5)*randn(5,10));

% change in covariance
x=[a;b];  
figure(2); clf
subplot 311
plot(x);
subplot 323
[idx,score_vector,score]=nt_split(x.^2); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('power')
subplot 324
[idx,score_vector,score]=nt_split(nt_xprod(x, 'lower')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('lower')
subplot 325
[idx,score_vector,score]=nt_split(nt_xprod(x, 'nodiag')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('nodiag')
subplot 326
[idx,score_vector,score]=nt_split(nt_normrow(nt_xprod(x,'lower'))); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('lower, normrow')

% change in covariance, change in power
x=[a;b;2*b];  
figure(3); clf
subplot 311
plot(x);
subplot 323
[idx,score_vector,score]=nt_split(x.^2); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('power')
subplot 324
[idx,score_vector,score]=nt_split(nt_xprod(x, 'lower')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('lower')
subplot 325
[idx,score_vector,score]=nt_split(nt_xprod(x, 'nodiag')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('nodiag')
subplot 326
[idx,score_vector,score]=nt_split(nt_normrow(nt_xprod(x,'lower'))); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('lower, normrow')

