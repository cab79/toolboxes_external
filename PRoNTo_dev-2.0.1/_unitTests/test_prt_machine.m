% simple test harness for prt_machine
%
% $Id$

%% test setup
featuresAsKernelMatrix=false;
useMultipleKernels=false;
useSynthData=true;         % use synthetic data
% root of PRT mat
%p_PRTroot='/Volumes/cs-research/intelsys/intelsys0/green/pattern/testdata/MoAEpilot/october5_2011';
p_PRTroot='/Volumes/cs-research/intelsys/intelsys0/green/pattern/testdata/MoAEpilot/GUITests';
fn_PRTtoUse='PRT.mat'; % which PRT mat we want to use

%% generate data
if useSynthData==true
    Ntr=96;                % number of training vectors
    D=64*64*64;                    % dimensionality of the feature space
    Nte=30;                 % number of testing vectors
    classOffset=1;       % higher value = easier problem
    
    if useMultipleKernels==true
        
        t_d.train={[rand(Ntr/2,D); rand(Ntr/2,D)+classOffset]; [randn(Ntr/2,D); randn(Ntr/2,D)+classOffset] };
        t_d.test={[rand(Nte/2,D); rand(Nte/2,D)+classOffset]; [randn(Nte/2,D); randn(Nte/2,D)+classOffset] };
    else
        t_d.train={[rand(Ntr/2,D); rand(Ntr/2,D)+classOffset]};
        t_d.test={[rand(Nte/2,D); rand(Nte/2,D)+classOffset]};
    end
    
    % can also label class1=2 class2=1 to be perverse
    t_d.tr_targets=([ones(Ntr/2,1)*1 ; ones(Ntr/2,1)*2]);
    te_targets=([ones(Nte/2,1)*1 ; ones(Nte/2,1)*2]);
else
    load(fullfile(p_PRTroot,fn_PRTtoUse));
    
    if length(PRT.group.subject.modality.design.conds)>2
        error('code not ready for >2 class testing')
    end
    nClasses=2;
    
    
    % no need to do explicit masking - already done in feature preparation
    % load mask
    %VMi=spm_vol(PRT.masks.fname(1:end-2));
    %VM=spm_read_vols(VMi);
    % find non-zero voxels (within the mask)
    % nzidx=find(VM>0); % no support for logical indexing in file_array :(
    %clear VM VMi;
    
    sz=size(PRT.fas.dat); % could also check size(PRT.fs.id_mat,1)
    T=sz(1); D=sz(2);
    
    %%% retrieve scan indices
    % make shortcuts
    c1ScansIdx=PRT.group.subject.modality.design.conds(1).scans;
    c1BlocksId=PRT.group.subject.modality.design.conds(1).blocks;
    c2ScansIdx=PRT.group.subject.modality.design.conds(2).scans;
    c2BlocksId=PRT.group.subject.modality.design.conds(2).blocks;
    % flag used scans (maybe useless)
    c1Scans=zeros(1,T);
    c2Scans=c1Scans;
    c1Scans(c1ScansIdx)=1;
    c2Scans(c2ScansIdx)=1;
    % recover block indices in terms of scan indices
    blocks=zeros(2,T);
    blocks(1,c1ScansIdx)=c1BlocksId;
    blocks(2,c2ScansIdx)=c2BlocksId;
    
    figure; subplot(211); stem(blocks(1,:),'fill'); axis tight;
    xlabel('scan index'); ylabel('block'); grid on;
    title(['c1 - ' PRT.group.subject.modality.design.conds(1).cond_name]);
    subplot(212); stem(blocks(2,:),'fill'); axis tight;
    xlabel('scan index'); ylabel('block'); grid on;
    title(['c2 - ' PRT.group.subject.modality.design.conds(2).cond_name]);

    
    %%% lazy feature generation - stupid block-mean extractor
    nBlocks=([numel(unique(c1BlocksId)) numel(unique(c2BlocksId))]);
    X=zeros(sum(nBlocks),D);
    exidx=1;
    for c=1:nClasses
        fprintf('Class %d block',c);
        for b=1:nBlocks(c)
            fprintf(' %d',b);
            scansIdx=find(blocks(c,:)==b);
            X(exidx,:)=mean(PRT.fas.dat(scansIdx,:),1);
            %X(:,exidx)=X(:,exidx)./std(PRT.file_arrays.Y(nzidx,scansIdx),[],2);
            exidx=exidx+1;
        end
        fprintf('%s\n','.');
    end
    figure; imagesc(X); xlabel('voxel'); ylabel('example');
    title('block means');
    
    %%% generate train/test labels for a 2-fold CV
    Ntr1=ceil(nBlocks(1)/2); Nte1=floor(nBlocks(1)/2);
    Ntr2=ceil(nBlocks(2)/2); Nte2=floor(nBlocks(2)/2);
    tridx=false(sum(nBlocks),1);
    tridx(1:Ntr1)=true; tridx((nBlocks(1)+1):(nBlocks(1)+Ntr2))=true;
    t_d.tr_targets=([ones(Ntr1,1)*1 ; ones(Ntr2,1)*2]);
    te_targets=([ones(Nte1,1)*1 ; ones(Nte2,1)*2]);
    % split into train and tes
    t_d.train={X(tridx,:)};
    t_d.test={X(~tridx,:)};
end

% compute kernel if needed
if featuresAsKernelMatrix==true
    if useMultipleKernels==true
        t_d.test={t_d.test{1}*t_d.train{1}';t_d.test{2}*t_d.train{2}'};
        t_d.train={t_d.train{1}*t_d.train{1}';t_d.train{2}*t_d.train{2}'};
    else
        t_d.test={t_d.test{1}*t_d.train{1}'};
        t_d.train={t_d.train{1}*t_d.train{1}'};
    end
    t_d.use_kernel=true;
else
    t_d.use_kernel=false;
end

t_d.pred_type='classification';

%% plot dataset
figure;
subplot(221); imagesc(t_d.train{1}); title('TR 1');
xlabel('feature'); ylabel('example');
subplot(223); imagesc(t_d.test{1}); title('TE 1');
if useMultipleKernels==true
    subplot(222); imagesc(t_d.train{2}); title('TR 2');
    subplot(224); imagesc(t_d.test{2}); title('TE 2');
end


%% prepare machine
%myMachine.function='prt_machine_svm_bin';
myMachine.function='prt_machine_RT_bin';

if ~isempty(strfind(myMachine.function,'svm_bin'))
    if featuresAsKernelMatrix==true
        myMachine.args='-s 0 -t 4';
    else
        myMachine.args='-s 0 -t 0';
    end
else
    myMachine.args=[601];
end

testCov=[];

%% potentially annoy the code
% 1: by removing libSVM, we might causes the code to point at the biostats
%toolbox's version of svmtrain instead
%rmpath('/Users/Richiardi/_skool/matlabTools/libsvm/'); 
% 2: by removing biostats toolbox, we're 
%rmpath('/Applications/_sciEng/MATLAB_R2011a.app/toolbox/bioinfo/biolearning');

%% run
tic
% OLD fomat
%output = prt_machine(train,test,testCov,myLabs,myMachine,featuresAsKernelMatrix);
output = prt_machine(t_d,myMachine);
toc

%% eval 
figure; plot(te_targets,'go','LineWidth',3,'MarkerSize',10); hold on;
plot(output.predictions,'k+','LineWidth',2,'MarkerSize',10);
plot(output.func_val,'b--','LineWidth',2); 
legend('true labels','hard prediction','soft prediction');
grid on
acc=sum(output.predictions==te_targets)/numel(te_targets)
title([myMachine.function ' - ' num2str(acc,'%2.2f')],'Interpreter','none');