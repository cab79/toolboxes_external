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

DSSonERP(data{2},[1 2 3],ones(1,size(data{2},3)),5,[],'on'); % maximise repeatability