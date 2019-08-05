//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data.
data{
  int<lower=1> n; // number of observations/samples...
  //int<lower=0, upper=2> genotype[n]; // Genotype
  vector[n] genotype;
  vector[n] y; // the expression data

  vector[n] measProp2; // the "measured" cell type proportions (i.e. the point estimate of these)
  vector[n] sd2; // the measurement error for each proportion estimate
}

// The parameters accepted by the model.
parameters{

  // intercept
  real beta0;

  // NOTE: Update to work for any number of cell types (beta will need to be a vector etc).
  real beta1;
  real beta2;
  real beta3;

  // Estimated precision in the underlying cell types (NOTE: At some point, make >=2 cell types script, this will need to be a vector).
  // Also, "precision" is used here as a throwback to implementing this originally in JAGS. This shold be re-build using variance/SD instead.
  real<lower=0> prec1;
  real<lower=0> prec2;

  // unknown/unmeasured *true* vector of proportions. These unmeasured values are treated as parameters.
  vector[n] trueProp2;

}

transformed parameters{

  // NOTE: update to arbitrary numbers of cell types
  real sigma1; // variance in cell type 1
  real sigma2;

  // we can compute betaNormal here as the sum of beta1 and 3. This needs to work for arbitrary numbers of cell types.
  real betaNormal;
  betaNormal = beta1 + beta3;

  // calcualte variance from precision (this is a throwback to when I implemented these in JAGS and it wanted the precision, this can be changed.)
  sigma1 = sqrt(1/prec1); // estimate of underlying variance in cell type 1
  sigma2 = sqrt(1/prec2); // estimate of underlying variance in cell type 2
}

// The model to be estimated.
model{

  // MODEL PRIORS //

  // weaky informative priors for these parameters
  beta0 ~ normal(0.0, 10);
  beta1 ~ normal(0.0, 10);
  beta2 ~ normal(0.0, 10);
  beta3 ~ normal(0.0, 10);
  prec1 ~ gamma(1.0/10, 1.0/10);
  prec2 ~ gamma(1.0/10, 1.0/10);


  // priors for each trueProp2 (this is vectorized)
  trueProp2 ~ normal(.5, 5); // weakly informative prior on each unobserved "true" proportion estimate (this vector syntax assigne this prior to each element of trueProp2)


  // MODEL LIKLIHOOD //
  measProp2 ~ normal(trueProp2, sd2); // error model. (vectorized)

  // This can't be vectorized because you apparently can't exponentiate a vector (i.e. trueProp^2 breaks vectorization).
  // Would it be possible to code this such that "trueProp" is just written as a massive series of strings? Maybe there's a way around this.
  // NOTE: write this for arbitrary numbers of cell types and probably swap precicion for variance (see PDF).
  for(i in 1:n)
  {
    y[i] ~ normal(beta0 + beta1*genotype[i] + beta2*trueProp2[i] + beta3*(trueProp2[i]*genotype[i]),
    sqrt(((1-2*trueProp2[i]+trueProp2[i]^2)/prec1) + (trueProp2[i]^2/prec2))
    );
  }
}

