

data {
  int<lower=0> N;
  vector[N] x;
  //array[N] real y;
  vector[N] ss; // sample size per survey
  vector[N]  npos;  // number positive for mutant
  int<lower=0> N_group;
  vector[N] int<lower=0, upper=N_group> group_id;
}
parameters {
  real alpha_0; // intercept
  real beta_0; // slope
  real<lower=0> sigma_g; // global sigma
  array[N_group] vector[2] random_group; // equvivalent to an N*2 matrix
  vector<lower=0>[2] sigs; // sds for random intercepts and slopes
  corr_matrix[2] R; // correlation matrix R
}
transformed parameters {
  cov_matrix[2] Sigma; // VCV matrix
  Sigma = quad_form_diag(R, sigs);
  array[N] real mu;
  array[N] real theta;

  // Fixed effects + Random effects
  for (i in 1 : N) {
    mu[i] = alpha_0 + beta_0 * x[i] + random_group[group_id[i], 1]
    + random_group[group_id[i], 2] * x[i];
    theta[i] = exp(mu[i]) / (1+ exp(mu[i]));
  }


}

model {
  alpha_0 ~ normal(0, 10);
  beta_0 ~ normal(0, 10);
  sigma_g ~ normal(0, 5);
  sigs ~ normal(0, 5);
  R ~ lkj_corr(2.0);
  random_group ~ multi_normal(rep_vector(0, 2), Sigma);

  npos ~ binomial(ss, theta);
}
