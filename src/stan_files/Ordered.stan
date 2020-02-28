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

  n = (N - 1.0) / 2.0;
  b = s2 .* N;
  nplus = sum(n);
}

parameters {
  real<lower=0> tau;
  positive_ordered[k] lambda_unconstrained;
}

transformed parameters {
  vector[k] sds;
  simplex[k] rho;
  vector[k] lambda;
  for (i in 1:k) lambda[i] = lambda_unconstrained[index_vector[i]];

  rho = lambda / sum(lambda);
  sds = 1.0 ./ sqrt(rho * tau * k);
}

model{
  target += -log(tau);
  lambda_unconstrained ~ gamma(alpha, 1);

  // adjust prior
  target += log(tgamma(nr_ordered + 1)); // for order constraints
  target += lgamma(alpha * (k - nr_equal)) - sum(lgamma(rep_vector(alpha, k - nr_equal))); // for equality constraints

  if (!(priors_only == 1)) {
      target += (
        ((k - sum(N))/2.0) * log(2*pi()) +
        dot_product(rep_vector(-0.50, k), log(N)) +
        nplus * log(tau*k) + dot_product(n, log(rho)) - k*tau*dot_product(b/2.0, rho)
      );
  }
}
