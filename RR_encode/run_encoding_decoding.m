% Run encoding to learn a model on training data
% Forming reconstructions from the encoding models and prior information
% with Gaussian mixture models and semantic gating
%
% ------------------------------------------------------------------------
% regularized linear regression toolbox, December 2014
% (c) Sanne Schoenmakers, Tom Heskes, Marcel van Gerven
% Donders Institute for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%
% Please refer to the following paper:
% Schoenmakers, S., Güçlü, U., Van_gerven, M., & Heskes, T. (2014). Gaussian 
% mixture models and semantic gating improve reconstructions from human 
% brain activity. Name: Frontiers in Computational Neuroscience, 8, 173.
%
% Kindly report any suggestions or corrections to
% s.schoenmakers@donders.ru.nl



clear all
restoredefaultpath
addpath(genpath('C:\Data\Matlab\RR_encode'))
dbstop if error
plot_on=1;

% analysis type
decode_type = 'uni'; %uni, multi, sem, unsup 

% set input directory
inputdir = 'C:\Data\temp\char_data\';

% load character-images that where shown in MRI scanner (trials x pixels)
load([inputdir 'Y_brains.mat']);

% display images
if 0
    figure
    for i=1:size(Y,1)
        image(reshape(Y(i,:),56,56));
        pause
    end
end

% load selection for the trainset and testset
load([inputdir 'train_testset.mat']);

% load the labels of the images that where shown in the MRI scanner (trials x labels)
load([inputdir 'labels']);

% load the voxels with the braindata for each image (trials x voxels)
load ([inputdir 'XS03_V1.mat']);

try
    % load prior data (images x pixels)
    load([inputdir 'prior_brains']);
    
    % construct prior, one prior containing all letter categories
    prior = [priorA;priorB;priorI;priorN;priorR;priorS];
    
catch
    % create priors from stimulus set (not ideal)
    ui = unique(L);
    for i = 1:length(ui)
        subY = Y(L==ui(i),:);
        prior(700*(i-1)+(1:700),:) = subY(randsample(size(subY,1),700, true),:);
%         av_img(i,:) = mean(Y(L==ui(i),:));
% 
%         av_img(i,:) = reshape(imgaussfilt(reshape(av_img(i,:),56,56),1),1,[]);
% 
%         maxI = prctile(av_img(i,:),75);
%         av_img(i,av_img(i,:)<maxI) = 0;
%         av_img(i,av_img(i,:)>0) = 1;
%         if 0
%             figure
%             imagesc(reshape(av_img(i,:),56,56));
%         end
    end
    %prior = reshape(repmat(reshape(av_img,1,[]),700,1),700*length(ui),[]);
end

% normalize prior
prior = zscore(double(prior));

%%%%%%%%%%%%%%%ENCODING%%%%%%%%%%%%%%%%

% start with seed in random number generator
rng(1);

% number of lambdas to do cross-val on
df_num = 100; 

% number of folds in the traindata
folds = 5; 

% z determines which method to use. Either 1 or 0.
z = 1; 

% use only traindata to train on
Lab = L(trainidx,1);
Im = Y(trainidx,:)';
Bold = X(trainidx,:)';

%run ridge regression
% x-fold cross-validated on folds of train data only. Test data used only
% for decoding.
[Beta_hat, Sigma_hat, X_mean, X_stand_de, Y_mean, Y_stand_de] = ...
    cross_validation(df_num, folds, Lab, Im, Bold, z);

switch decode_type
    
    %decode_type
    %disp(['running decoding type: ' decode_type])
    
    case 'uni'
        %%%%%%%%%%%% UNIMODAL DECODING %%%%%%%%%%%%%%%%%%%

        % calculate the covariance of the prior
        % Prior is a multivariate Gaussian with zero mean (z-scored above), so only need to consider covariance. 
        R = cov(prior);
        figure; imagesc(R);colormap('hot');

        % run unimodal decoding
        % train_idx fed into this function to calculate mean/std of X/Y in
        % training set in order to normalise test set.
        [images, recons] = decode_unimod(X, Y, trainidx, testidx, Beta_hat', Sigma_hat', R);

        % plot images and their reconstructions
        if plot_on
            for j=1:size(recons,1) 
                subplot(1,2,1); imagesc(reshape(images(j,:),[56 56]));    
                axis off; 
                title(j);
                subplot(1,2,2); imagesc(reshape(recons(j,:),[56 56]));     
                axis off; 
                colormap('hot');
                pause(1); 
            end
        end

        %%%%%%%%%%%% MULTIMODAL DECODING %%%%%%%%%%%%%%%%%%%
    case 'multi'
        % create prior with multiple priorsets
        multipriordata = zeros(700,3136,6);
        for i=1:6
          multipriordata(:,:,i) = prior(700*(i-1)+(1:700),:);
        end

        %set temperature for mixture in result, T=1 results in an equal mixture of
        %all components, T approaches 0 results in zooming in on the most likely
        %category
        T = [0.001 1];

        %decode with multimodal prior
        [images, recons, w, n, mu] = decode_mixtmod(X, Y, trainidx, testidx, Beta_hat', Sigma_hat', multipriordata,T);

        % plot images and their reconstructions
        if plot_on
            for Temp = 1: size(T,2)
                for j=1:size(recons{1},1) 
                    subplot(1,2,1); imagesc(reshape(images(j,:),[56 56]));    
                    axis off; 
                    title(j);
                    subplot(1,2,2); imagesc(reshape(recons{Temp}(j,:),[56 56]));     
                    axis off; 
                    colormap('hot');
                    pause(1); 
                end
            end
        end

        %%%%%%%%%%%% MULTIMODAL DECODING WITH SEMANTIC GATING %%%%%%%%%%%%%%%%%%%
    case 'sem'
        % load the voxels with the braindata for each image
        load ([inputdir 'XS03_V2.mat']);

        %set temperature for mixture in result, T=1 results in an equal mixture of
        %all components, T approaches 0 results in zooming in on the most likely
        %category
        T = [0.001 1];

        % start with seed in random number generator
        rng(1);

        %add functions to matlabpath
        addpath ./minFunc/
        addpath ./glmnet_matlab/

        %chosse number of folds
        folds = 10;

        %set sparsity, lasso = sparse = 1, ridge = 0 
        is_sparse = 1;

        %Run multinomial regression to get probability of each category 
        [alpha_hat_V2, Beta_hat_V2] = cross_validation2(is_sparse, folds, X_V2(trainidx,:)', L(trainidx,:)');
        [Pmlr, y_hat_V2] = validation(Beta_hat_V2, X_V2(testidx,:)');

        %Decode with multimodal prior and semantic gating
        [images, recons, w, n, mu] = decode_mixtmod_mlr(X, Y, trainidx, testidx, Beta_hat', Sigma_hat', multipriordata, T, Pmlr);

        % plot images and their reconstructions
        if plot_on
            for Temp = 1: size(T,2)
                for j=1:size(recons{1},1) 
                    subplot(1,2,1); imagesc(reshape(images(j,:),[56 56]));    
                    axis off; 
                    title(j);
                    subplot(1,2,2); imagesc(reshape(recons{Temp}(j,:),[56 56]));     
                    axis off; 
                    colormap('hot');
                    pause(1); 
                end
            end
        end

    case 'unsup'
        %%%%%%%%%%%% UNSUPERVISED MULTIMODAL DECODING %%%%%%%%%%%%%%%%%%%

        %set temperature for mixture in result, T=1 results in an equal mixture of
        %all components, T approaches 0 results in zooming in on the most likely
        %category
        T = [0.001 1];

        %set number of clusters
        clust = 20;

        %we start with a random seed and run it several times, each seed gives a different result
        rng('shuffle')

        %discart empty pixels in the prior
        sel = any(prior);
        pixnum = sum(sel);

        %Calculate which instances go in which cluster
        [IDX,mu_prior] = kmeans(prior(:,sel),clust);

        %decode with an unsupervised multimodal prior
        [images, recons, w, n, mu] = decode_mixtmod_manyK(X, Y, trainidx, testidx, Beta_hat', Sigma_hat', mu_prior, IDX, prior, sel, T);

        % plot images and their reconstructions
        if plot_on 
            for Temp = 1: size(T,2)
                for j=1:size(recons{1},1) 
                    subplot(1,2,1); imagesc(reshape(images(j,:),[56 56]));    
                    axis off; 
                    title(j);
                    subplot(1,2,2); imagesc(reshape(recons{Temp}(j,:),[56 56]));     
                    axis off; 
                    colormap('hot');
                    pause(1); 
                end
            end
        end

end