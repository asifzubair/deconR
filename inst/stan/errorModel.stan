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
  // number of cell types
  int<lower=0> numCellTypes;
  // number of observations/samples...
  int<lower=1> numSamples;
  //int<lower=0, upper=2> genotype[n]; // Genotype
  vector[numSamples] genotype;
  // the expression data
  vector[numSamples] y;

  // the "measured" cell type proportions (i.e. the point estimate of these)
  vector[numSamples] measProp;
  // the measurement error for each proportion estimate
  vector[numSamples] sd2;
}

// The parameters accepted by the model.
parameters{

  // intercept
  real beta0;

  // NOTE: Update to work for any number of cell types (beta will need to be a vector etc).
  vector[numCellTypes] beta;

  // Estimated precision in the underlying cell types
  // (NOTE: At some point, make >=2 cell types script, this will need to be a vector).
  // Also, "precision" is used here as a throwback to implementing this originally in JAGS.
  // This shold be re-build using variance/SD instead.
  vector<lower=0>[numCellTypes] prec;

  // unknown/unmeasured *true* vector of proportions. These unmeasured values are treated as parameters.
  vector[numSamples] trueProp;

}

transformed parameters{

  // NOTE: update to arbitrary numbers of cell types
  // variance in cell type 1
  vector<lower=0>[numCellTypes] sigma;
  real sigma2;

  // we can compute betaNormal here as the sum of beta1 and 3.
  // This needs to work for arbitrary numbers of cell types.
  real betaNormal;
  betaNormal = beta[1] + beta[3];

  // calcualte variance from precision
  // (this is a throwback to when I implemented these in JAGS and it wanted the precision, this can be changed.)
  sigma[1] = sqrt(1/prec[1]); // estimate of underlying variance in cell type 1
  sigma[2] = sqrt(1/prec[2]); // estimate of underlying variance in cell type 2
}

// The model to be estimated.
model{

  // MODEL PRIORS //

  // weaky informative priors for these parameters
  beta0 ~ normal(0.0, 10);
  beta[1] ~ normal(0.0, 10);
  beta[2] ~ normal(0.0, 10);
  beta[3] ~ normal(0.0, 10);
  prec[1] ~ gamma(1.0/10, 1.0/10);
  prec[2] ~ gamma(1.0/10, 1.0/10);


  // priors for each trueProp2 (this is vectorized)
  // weakly informative prior on each unobserved "true" proportion estimate
  // (this vector syntax assigne this prior to each element of trueProp2)
  trueProp2 ~ normal(.5, 5);


  // MODEL LIKELIHOOD
  // error model. (vectorized)
  measProp ~ normal(trueProp2, sd2);

  // This can't be vectorized because you apparently can't exponentiate a vector (i.e. trueProp^2 breaks vectorization).
  // Would it be possible to code this such that "trueProp" is just written as a massive series of strings?
  // Maybe there's a way around this.
  // NOTE: write this for arbitrary numbers of cell types and probably swap precicion for variance (see PDF).
  for(i in 1:numSamples)
  {
    y[i] ~ normal(beta0 + beta[1]*genotype[i] + beta[2]*trueProp2[i] + beta[3]*(trueProp2[i]*genotype[i]),
    sqrt(((1-2*trueProp2[i]+trueProp2[i]^2)/prec[1]) + (trueProp2[i]^2/prec[2]))
    );
  }
}

