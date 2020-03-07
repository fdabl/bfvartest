// This implements the only-order-constraint model, accommodating hypotheses such as '1<2<3<4'
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
  // real ll_const; // constant from log likelihood

  n = (N - 1.0) / 2.0;
  b = s2 .* N;
  nplus = sum(n);
  // ll_const = -0.5 * sum(log(N)) + (k - sum(N)) / 2.0 * log(2 * pi());
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

  // adjust prior
  target += lgamma(k + 1);

  if (!(priors_only == 1)) {
    // target += ll_const;
    target +=  dot_product(n, log(prec));
    target += -0.5 * dot_product(prec, b);
  }
}
