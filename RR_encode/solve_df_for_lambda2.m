function lamb = solve_df_for_lambda2(d, num_of_df)
%SOLVE_DF_FOR_LAMBDA solves the degrees of freedom of the ridge regression
%fit for the ridge parameter (lambda)
%
%   INPUT
%
%   d is the singular values of the input matrix
%   num_of_df is the number of degrees of freedom to solve for lambda
%
%   OUTPUT
%
%   lamb is the estimates of lambda
%

guess       = @(d, df, lamb, N) max(1 / mean(d .^ -2) * (N - df) / df, lamb);
f           = @(d, df, lamb)  sum(d .^ 2 ./ (d .^ 2 + lamb)) - df;
f_prime     = @(d, lamb) -sum(d .^ 2 ./ (d .^ 2 + lamb) .^ 2);
iterate     = @(d, df, lamb) max(0, lamb - f(d, df, lamb) / f_prime(d, lamb));
length_of_d = length(d);
df          = linspace(1, length_of_d, num_of_df);

disp(['df: ' num2str(num_of_df)]);

lamb(num_of_df) = guess(d, df(num_of_df), 0 , length_of_d);

while abs(f(d, df(num_of_df), lamb(num_of_df))) > 1e-12
    
    lamb(num_of_df) = iterate(d, df(num_of_df), lamb(num_of_df));
    
end

for i = num_of_df - 1 : -1 : 1
    
    disp(['df: ' num2str(i)]);
    
    lamb(i) = guess(d, df(i), lamb(i + 1), length_of_d);

    while abs(f(d, df(i), lamb(i))) > 1e-12
        
        lamb(i) = iterate(d, df(i), lamb(i));
        
    end
    
end

end

