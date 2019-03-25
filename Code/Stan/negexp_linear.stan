// Thanks to Ben Goodrich: https://discourse.mc-stan.org/t/help-translating-model-from-bugs/8198
data {
  int<lower = 0> N_obs;   // observations on x (climate)
  vector[N_obs] x_obs;
  int<lower = 0> N_miss; // missing observations on x to estimate
  vector[N_obs + N_miss] log_y;
  vector[N_obs + N_miss] age;

  // I'm not really sure what these are
  int<lower = 0> M; // number of trees?
  int<lower = 1, upper = M> tree_ID[N];
}
transformed data {
  int N = N_obs + N_miss;
}
parameters {
  vector[N_miss] x_miss;
  vector[N] eta_raw;

  vector[M] alpha_0_raw;
  vector<upper = 0>[M] alpha_1;
  vector<lower = 0>[M] sigma_y;

  // need good priors for all of these
  real mu_x;
  real<lower = 0> sigma_x;
  real<lower = 0> sigma_eta;
  real beta;
  real mu_a0;
  real<upper = 0> mu_a1;
  real<lower = 0> sd_a0;
  real<lower = 0> sd_a1;
}
model {
  vector[N] x = append_row(x_miss, x_obs);
  vector[N] eta = beta * x + sigma_eta * eta_raw;
  vector[N] alpha0 = mu_a0 + sd_a0 * alpha0_raw; 
 
  log_y ~ normal(alpha_0[tree_ID] + alpha_1[tree_ID] * age + eta, sigma_y);
  x ~ normal(mu_x, sigma_x);
  eta_raw ~ std_normal(); // implies eta ~ normal(beta * x, sigma_eta);
  
  alpha0_raw ~ std_normal(); // implies alpha0 ~ normal(mu_a0, sd_a0)
  alpha1 ~ normal(mu_a1, sd_a1);
  target += -M * normal_lccdf(0 | mu_a1, sd_a1); // accounts for truncation
  sd_y ~ student_t(3, 0, 0.02); // do not need to account for truncation
  
  // put good priors here on all those scalars in the parameters block
  sigma_x ~ cauchy(0, 2.5);
  sigma_eta ~ cauchy(0, 2.5);
  beta ~ normal(0, 3);
  mu_a0 ~ normal(6, 3); // weakly informative prior with mean log(y) as mean
  mu_a1 ~ normal(0, 3);
  sd_a0 ~ cauchy(0, 2.5);
  sd_a1 ~ cauchy(0, 2.5);
}