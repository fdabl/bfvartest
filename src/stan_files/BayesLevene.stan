functions {

  vector order_parameters(vector rho, vector beta, int k, int nr_rel, int[] relations, int[,] constraint_mat) {
    vector[k] sorted_rho = rho;

    int is_ordered = 1;
    int is_unconstrained = 1;
    for (i in 1:nr_rel) {
      if (relations[i] != 2) {
        is_ordered = 0;
      }
    }

    for (i in 1:nr_rel) {
      if (relations[i] != 3) {
        is_unconstrained = 0;
      }
    }

    if (!(is_unconstrained == 1)) {
      {
        int i;
        int bi;
        int hi;
        int lo;
        int el;
        int cur;
        int rel;
        int prev;
        vector[k] ubound;
        vector[k] lbound;
        int nr_equals;
        int nr_equalities;
        int relation_between;
        int cur_constraint;
        int nr_unconstraint;
        bi = 1;
        el = 1;
        cur_constraint = 1;

        // loop over each element, and adjust it according to the hypothesis
        while (el <= k) {

          // find the upper and lower constraints of the current element
          lo = constraint_mat[cur_constraint, 1];
          hi = constraint_mat[cur_constraint, 2];

          if (lo == -99) {
            lbound = rep_vector(0, k);
          } else {
            lbound = rep_vector(0, k);
            lbound[lo:k] = sorted_rho[lo:k];
          }

          if (hi == -99) {
            ubound = rep_vector(1, k);
          } else {
            ubound = rep_vector(1, k);
            ubound[1:hi] = sorted_rho[1:hi];
          }

          // the element has NO constraints!
          if (min(ubound) == 1.0 && max(lbound) == 0.0) {

            // if the relationship is NOT unconstrained (e.g., 1,2) then do update!
            if (relations[el] != 3) {
              sorted_rho[el] = beta[bi] * (min(ubound) - max(lbound)) + max(lbound);
            }

          } else {
            sorted_rho[el] = beta[bi] * (min(ubound) - max(lbound)) + max(lbound);
          }

          // check whether elements are equal, and if so, how many
          i = 0;
          nr_equalities = 0;
          while ((el + i < k) && relations[el + i] == 1) {
            nr_equalities = nr_equalities + 1;
            i = i + 1;
          }

          // bulk update the contraints for elements that are equal
          for (j in 1:nr_equalities) {
            sorted_rho[el + j] = beta[bi] * (min(ubound) - max(lbound)) + max(lbound);
          }

          bi = bi + 1;
          el = el + nr_equalities; // fast forward to elements that are not equal

          // if next one is an order constraint, update the upper and lower constraints
          if ((el + 1) <= k && relations[el] == 2) {
              cur_constraint = cur_constraint + 1;
          }

          el = el + 1;
        }
      }
    }

    sorted_rho = sorted_rho / sum(sorted_rho);
    return sorted_rho;
  }
}

data {
  int k;
  real alpha;
  int nr_rel;
  int nr_free;
  int nr_equal;
  int nr_ordered;
  int index_vector[k];
  int relations[nr_rel];
  int constraint_mat[nr_ordered + 1, 2];
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
  vector<lower=0, upper=1>[k] beta;
  simplex[k] rho_unconstrained;
}

transformed parameters {
  simplex[k] rho;

  for (i in 1:k) {
    rho[i] = rho_unconstrained[index_vector[i]];
  }

  rho = rho / sum(rho);

  // "1,2>3>4"
  // rho[1] = beta[1] * (1 - max(rho[3:4])) + max(rho[3:4]);
  // rho[2] = beta[2] * (1 - max(rho[3:4])) + max(rho[3:4]);
  // rho[3] = beta[3] * (min(rho[1:2]) - rho[4]) + rho[4];
  // rho[4] = beta[4] * (min(rho[1:3]) - 0) + 0;

  rho = order_parameters(rho, beta, k, nr_rel, relations, constraint_mat);
}

model{
  target += -log(tau);
  target += dirichlet_lpdf(rho_unconstrained | rep_vector(alpha, k));

  if (!(priors_only == 1)) {
      target += (
        ((k - sum(N))/2) * log(2*pi()) +
        dot_product(rep_vector(-0.50, k), log(N)) +
        nplus * log(tau*k) + dot_product(n, log(rho)) - k*tau*dot_product(b/2, rho)
      );
  }
}
