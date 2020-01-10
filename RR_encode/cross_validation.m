function [Be_hat, Sig_hat, X_mean, X_stand_de, Y_mean, Y_stand_de] = ...
    cross_validation(df_num, K, N, X, Y, z)
%CROSS_VALIDATION Estimates a ridge regression model using cross-validation
%
%   INPUT
%
%   df_num is the number of ridge parameters (lambdas) to test
%   K is the number of cross-validation folds
%   N is the number of trials
%   X is the training set (images) - X should be p * N where p is number of
%   pixels and N is number of images
%   Y is the training set (voxels) - Y should be v * N where v is number of
%   voxels and N is number of images
%   z is a boolean to indicate whether to rotate the data (z == true) to get
%   rid of zero eigenvalues or not (z == false)
%
%   OUTPUT
%
%   Be_hat is the estimated regression parameters
%   Sig_hat is the estimated noise covariance matrix
%   X_mean, X_stand_de, Y_mean, Y_stand_de are the mean and standard
%   deviations of the input and output
%

[X, X_mean, X_stand_de] = zscore(X, [], 2);
[Y, Y_mean, Y_stand_de] = zscore(Y, [], 2);

x_val_in  = crossvalind('Kfold', N, K);
[U, D, V] = svd(X, 'econ');
d         = diag(D);
lamb_hat  = solve_df_for_lambda2(d, df_num); % a
 %lamb_hat  = logspace(5, -5, df_num); % b


for i = K : -1 : 1
    
    
    val_in    = x_val_in == i;
    es_in     = ~val_in;
    
    if z == 1
    
    [U_es, D_es, V_es] = svd(X(:, es_in), 'econ');
    d_es         = diag(D_es);
    d_es_squared = d_es .^ 2;
    
    end
    
    for j = df_num : -1 : 1
        disp(['k: ' num2str(i) ', j: ' num2str(j)]);
        
        if z == 1
        
        re_var{i}(j, :, :) = (Y(:, val_in) - Y(:, es_in) * V_es * ...
            diag(d_es ./ (d_es_squared + lamb_hat(j))) * U_es' * X(:, ...
            val_in))';
        
        else
        
        re_var{i}(j, :, :) = (Y(:, val_in) - Y(:, es_in) / (X(:, ...
            es_in)' * X(:, es_in) + lamb_hat(j) * eye(sum(es_in))) * ...
            X(:, es_in)' * X(:, val_in))';
        
        end
        
    end
    
end

re_var                      = squeeze(var(cell2mat(re_var), [], 2));
[re_var_min, re_var_min_in] = min(re_var);
Y_var                       = var(Y, [], 2) - 0.01;

if z == 1

d_squared = d .^ 2;

end

X_size  = size(X);

for i = size(Y, 1) : -1 : 1
    
    disp(['Calculating betas: ' num2str(i)]);
    
    if re_var_min(i) < Y_var(i)
        
        if z == 1
        
        Be_hat(i, :)  = Y(i, :) * V * diag(d ./ (d_squared + ...
            lamb_hat(re_var_min_in(i)))) * U';
        
        else
        
        Be_hat(i, :)  = Y(i, :) / (X' * X + ...
            lamb_hat(re_var_min_in(i)) * eye(X_size(2))) * X';
        
        end
        
        Sig_hat(i, i) = re_var_min(i);
    
    else
        
        Be_hat(i, :)  = zeros(1, X_size(1));
        Sig_hat(i, i) = 1;
        
    end
    
end

end

