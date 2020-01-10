function [alpha_hat, Beta_hat] = cross_validation(is_sparse, K, X, y)
%CROSS_VALIDATION Estimates a multinomial logistic regression model  using
%cross-validation
%
%   INPUT
%
%   is_sparse is a Boolean for indicating whether to use lasso (is_sparse
%   == true) or ridge (is_sparse == false)
%   K is the number of folds
%   X is a p by N input matrix
%   y is a 1 by N label vector (with integers elements running from 1 to
%   number of categories)
%
%   OUTPUT
%
%   alpha_hat is the optimal regularization parameter
%   Beta_hat is the model estimated using alpha_hat
%

if ~ismatrix(X),                    error('~ismatrix(X)');                    end
if ~isvector(y),                    error('~isvector(y)');                    end
if size(X, 2) ~= size(y, 2),        error('size(X, 2) ~= size(y, 2)');        end
if ~isequal(unique(y), 1 : max(y)), error('~isequal(unique(y), 1 : max(y))'); end

cross_validation_indices = crossvalind('Kfold', y, K);

if is_sparse == 1
    
    fit = cvglmnet(X', y', 'multinomial', [], [], [], cross_validation_indices);
    
    cvglmnetPlot(fit);
    
    alpha_hat = fit.lambda_1se;
    
    for i = length(fit.glmnet_fit.beta) : -1 : 1
        
        fit.glmnet_fit.beta{i} = fit.glmnet_fit.beta{i}(:, alpha_hat == fit.glmnet_fit.lambda);
        
    end
    
    fit.glmnet_fit.a0 = fit.glmnet_fit.a0(:, alpha_hat == fit.glmnet_fit.lambda);
    
    Beta_hat = [cell2mat(fit.glmnet_fit.beta)' fit.glmnet_fit.a0];
    
    return;
    
end

options.alpha = 0;
fit           = glmnet(X', y', 'multinomial', options);
alpha         = fit.lambda;

for i = K : -1 : 1
    
    validation_indices = cross_validation_indices == i;
    estimation_indices = ~validation_indices;
    
    for j = length(alpha) : -1 : 1
        
        [~, y_hat]                    = validation(estimation(alpha(j), X(:, estimation_indices), y(estimation_indices)), X(:, validation_indices));
        classification_accuracy(i, j) = mean(y_hat == y(validation_indices));
        
    end
    
end

[~, index_of_alpha_hat] = max(mean(classification_accuracy));

alpha_hat = alpha(index_of_alpha_hat);
Beta_hat  = estimation(alpha_hat, X, y);

end

