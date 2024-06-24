data {
    int<lower=0> J; // number of schools
    array[J] real y; // estimated treatment effects
    array[J] real<lower=0> sigma; // s.e. of effect estimates
}
parameters {
    real mu;
    real<lower=0> tau;
    array[J] real eta;
}
transformed parameters {
    array[J] real theta;
    for (j in 1:J)
        theta[j] = mu + tau * eta[j];
}
model {
    eta ~ normal(0, 1);
    y ~ normal(theta, sigma);
}