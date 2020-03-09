// This implements the mixed constrained model.
// For example, if '1,2<3=4,5<6', then we sample from the reduced model '1,2,3,4,5'.
// If, on the other hand, the model does not have equalities, we sample from the full, unconstrained model.
data {
  int k;
  real alpha;
  int nr_equal;
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
  real<lower=0> lambda_unconstrained[k - nr_equal];
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

  if (!(priors_only == 1)) {
    // target += ll_const;
    target +=  dot_product(n, log(prec));
    target += -0.5 * dot_product(prec, b);
  }
}
