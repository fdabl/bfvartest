// This implements the only-order-constraint model, such as for example '1<2<3<4'
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

  n = (N - 1.0) / 2.0;
  b = s2 .* N;
  nplus = sum(n);
}

parameters {
  real<lower=0> tau;
  positive_ordered[k] lambda_unconstrained;
}

transformed parameters {
  vector[k] lambda;
  simplex[k] rho;
  for (i in 1:k) lambda[i] = lambda_unconstrained[index_vector[i]];

  rho = lambda / sum(lambda);
}

model{
  target += -log(tau);
  // lambda_unconstrained ~ gamma(alpha, 1);
  target += gamma_lpdf(lambda_unconstrained | alpha, 1);

  // // add normalizing constant of the reduced dirichlet prior
  // target += sum(lgamma(rep_vector(alpha, k - nr_equal))) - log(alpha * (k - nr_equal));

  if (!(priors_only == 1)) {
      target += (
        ((k - sum(N))/2) * log(2*pi()) +
        dot_product(rep_vector(-0.50, k), log(N)) +
        nplus * log(tau*k) + dot_product(n, log(rho)) - k*tau*dot_product(b/2, rho)
      );
  }
}
