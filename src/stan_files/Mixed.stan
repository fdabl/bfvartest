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

  n = (N - 1.0) / 2.0;
  b = s2 .* N;
  nplus = sum(n);
}

parameters {
  real<lower=0> tau;
  simplex[k] rho_unconstrained;
}

transformed parameters {
  simplex[k] rho;
  for (i in 1:k) rho[i] = rho_unconstrained[index_vector[i]];
  rho = rho / sum(rho);
}

model{
  target += -log(tau);
  rho_unconstrained ~ dirichlet(rep_vector(alpha, k));

  // adjust prior for equality constraints
  target += lgamma(alpha * (k - nr_equal)) - sum(lgamma(rep_vector(alpha, k - nr_equal)));

  if (!(priors_only == 1)) {
      target += (
        ((k - sum(N))/2) * log(2*pi()) +
        dot_product(rep_vector(-0.50, k), log(N)) +
        nplus * log(tau*k) + dot_product(n, log(rho)) - k*tau*dot_product(b/2, rho)
      );
  }
}
