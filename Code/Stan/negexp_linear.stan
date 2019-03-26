// Thanks to Ben Goodrich: https://discourse.mc-stan.org/t/help-translating-model-from-bugs/8198
data {
  int<lower = 0> N_obs;   // observations on x (climate)
  vector[N_obs] x_obs;
  int<lower = 0> N_miss; // missing observations on x to estimate
  vector[N_obs + N_miss] log_y;
  vector[N_obs + N_miss] age;

  // I'm not really sure what these are
  int<lower = 0> M; // number of trees
  int<lower = 0> N;
  int<lower = 1, upper = M> tree_id[N];
  int<lower = 1, upper = N> year[N];
}
transformed data {
  // int N = N_obs + N_miss;
}
parameters {
  vector[N_miss] x_miss;
  vector[N] eta_raw;

  vector[N] alpha0_raw;
  vector<upper = 0>[N] alpha1;
  vector<lower = 0>[N] sd_y;
  
 //  vector[M] alpha0_raw;
 //  vector<upper = 0>[M] alpha1;
 // // real<upper = 0> alpha1;
 //  vector<lower = 0>[M] sd_y;

  // need good priors for all of these
  real mu_x;
  real<lower = 0> sd_x;
  real<lower = 0> sd_eta;
  real beta0;
  real mu_a0;
  real<upper = 0> mu_a1;
  real<lower = 0> sd_a0;
  real<lower = 0> sd_a1;
}
model {
  vector[N] x = append_row(x_miss, x_obs)[year];
  vector[N] eta = beta0 * x + sd_eta * eta_raw;
  vector[N] alpha0 = mu_a0 + sd_a0 * alpha0_raw; 
 
 log_y ~ normal(alpha0[tree_id] + alpha1[tree_id] .* age + eta, sd_y[tree_id]);
 // for(n in 1:N) {
 // log_y[n] ~ normal(alpha0[tree_id[n]] + alpha1[tree_id[n]] * age[n] + eta[n], sd_y[tree_id[n]]);
 // }
  // log_y ~ normal(alpha0[tree_id] + alpha1 * age + eta, sd_y);
  x ~ normal(mu_x, sd_x);
  eta_raw ~ std_normal(); // implies eta ~ normal(beta * x, sigma_eta);
  
  alpha0_raw ~ std_normal(); // implies alpha0 ~ normal(mu_a0, sd_a0)
 alpha1 ~ normal(mu_a1, sd_a1);
 target += -M * normal_lccdf(0 | mu_a1, sd_a1); // accounts for truncation
 // sd_y ~ student_t(3, 0, 0.02); // do not need to account for truncation
  
  // put good priors here on all those scalars in the parameters block
  sd_x ~ cauchy(0, 2.5);
  sd_eta ~ cauchy(0, 2.5);
  beta0 ~ normal(0, 3);
  mu_a0 ~ normal(6, 3); // weakly informative prior with mean log(y) as mean
  mu_a1 ~ normal(0, 3);
  sd_a0 ~ cauchy(0, 2.5);
  sd_a1 ~ cauchy(0, 2.5);
}
