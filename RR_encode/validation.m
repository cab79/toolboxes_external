function [Pr, y_hat] = validation(Beta_hat, X)
%VALIDATION validates the MLR model estimated using cross_validation2.m
%
%   INPUT
%
%   Beta_hat is the output of the cross_validation2.m
%   X is the test set (n dims x n trials)
%
%   OUTPUT
%
%   Pr is the probability distribution across C categories for each of the
%   n trials, i.e. a C x n matrix
%   y_hat is the predicted labels for each of the n trials
%

if ~ismatrix(Beta_hat), error('~ismatrix(Beta_hat)'); end

M_plus_1 = size(X, 1) + 1;

if ~isequal(size(Beta_hat, 2), M_plus_1), error('~isequal(size(Beta_hat, 2), size(X, 1) + 1)'); end

X(M_plus_1, :) = 1;

Pr          = Beta_hat * X;
Pr          = exp(bsxfun(@minus, Pr, max(Pr)));
Pr          = bsxfun(@rdivide, Pr, sum(Pr));
[~, y_hat]  = max(Pr);

end

