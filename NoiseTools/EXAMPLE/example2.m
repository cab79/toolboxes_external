% Same as example1, but the data now include multiple conditions.
% We look for the linear combination that maximizes repeatability jointly
% for all conditions.  Data are in a cell array of matrices of dimensions 
% time*channels*trials
%
% Uses nt_dss0().

clear;
disp(mfilename);
help(mfilename)

% create synthetic data
nsamples=100*3;
nchans=30;
ntrials=100;
noise_dim=20; % dimensionality of noise
freqs=[1 2];
mix1=randn(1,nchans);
mix2=randn(noise_dim,nchans);
for iCondition=1:2
    source{iCondition}=[zeros(nsamples/3,1);sin(2*pi*freqs(iCondition)*(1:nsamples/3)/(nsamples/3))';zeros(nsamples/3,1)]; 
    s=source{iCondition}*mix1;
    s=repmat(s,[1,1,100]); % evoked
    SNR=0.1;
    noise=nt_mmat(randn(nsamples,noise_dim,ntrials), mix2);
    data{iCondition}=noise/rms(noise(:))+SNR*s/nt_rms(s(:));
end

% apply DSS to clean them
nKEEP = []; % leave empty to only remove default noisest components (see nt_dss0 defaults)
c0=zeros(nchans); c1=zeros(nchans);
for iCondition=1:2
    c0=c0+nt_cov(data{iCondition});
    c1=c1+nt_cov(mean(data{iCondition},3));
end
[todss,pwr0,pwr1]=nt_dss0(c0,c1,nKEEP);
for iCondition=1:2
    z{iCondition}=nt_mmat(data{iCondition},todss);
    newdata{iCondition} = nt_mmat(z{iCondition},todss');
end

% plot results
figure(1); clf
subplot 241; 
plot(source{1}); title('source, condition 1'); 
subplot 242;
plot(mean(data{1},3)); title('data');
subplot 243;
nt_bsplot(z{1}(:,1,:)); title('recovered'); 
subplot 244;
plot(mean(newdata{1},3)); title('reconstructed'); 
subplot 245; 
plot(source{2}); title('source, condition 2'); 
subplot 246;
plot(mean(data{2},3)); title('data');
subplot 247;
nt_bsplot(z{2}(:,1,:)); title('recovered'); 
subplot 248;
plot(mean(newdata{2},3)); title('reconstructed'); 

