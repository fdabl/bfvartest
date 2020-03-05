// This implements the only-order-constraint model, accommodating hypotheses such as '1<2<3<4'
functions {
  real log_prior_order_and_equality_constraints(int nr_ordered, int nr_equal, int k, real alpha) {

    real ans;
    ans  = lgamma(nr_ordered + 1);
    ans += lgamma(alpha * (k - nr_equal)) - (k - nr_equal) * lgamma(alpha);

    return ans;
  }
}
data {
  int k;
  real alpha;
  int nr_equal;
  int nr_ordered;
  int index_vector[k];
  vector<lower = 0>[k] s2;
  vector<lower = 0>[k] N;
  int priors_only;
}

transformed data {
  vector<lower = 0>[k] n;
  vector<lower = 0>[k] b;
  real nplus;
  real ll_const; // constant from log likelihood
  real lp_const; // constant from log prior

  n = (N - 1.0) / 2.0;
  b = s2 .* N;
  nplus = sum(n);
  ll_const = -0.5 * sum(log(N)) + (k - sum(N)) / 2.0 * log(2 * pi());
  lp_const = log_prior_order_and_equality_constraints(nr_ordered, nr_equal, k, alpha);
}

parameters {
  real<lower=0> tau;
  positive_ordered[k] lambda_unconstrained;
}

transformed parameters {
  vector[k] sds;
  simplex[k] rho;
  vector[k] lambda;
  vector[k] prec;
  for (i in 1:k) lambda[i] = lambda_unconstrained[index_vector[i]];

  rho = lambda / sum(lambda);
  prec= rho * tau * k;
  sds = 1.0 ./ sqrt(rho * tau * k);
}

model{
  target += -log(tau);
  lambda_unconstrained ~ gamma_lpdf(alpha, 1);

  target += lgamma(k + 1);

  // target += -sum(lambda_unconstrained);

  // adjust prior
  // target += lp_const;
  // target += log(tgamma(nr_ordered + 1)); // for order constraints
  // target += lgamma(alpha * (k - nr_equal)) - sum(lgamma(rep_vector(alpha, k - nr_equal))); // for equality constraints

  if (!(priors_only == 1)) {
    // target += ll_const;
    target +=  dot_product(n, log(prec));
    target += -0.5 * dot_product(prec, b);
  }
}
