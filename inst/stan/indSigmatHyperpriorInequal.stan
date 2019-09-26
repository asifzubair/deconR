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


// Following the logic from here:
// https://discourse.mc-stan.org/t/how-to-add-inequality-constraint-on-sum-of-parameters/1414 
// and here:
// https://discourse.mc-stan.org/t/how-to-add-inequality-constraint-on-the-model/4835/2 
// which basically recommends to downscale a simplex

// The input data
data
{
  int<lower=0> numGenes; // The number of "genes" (1000)
  int<lower=0> numCellTypes; // number of cell types (2)
  vector[numGenes] exprMixVec; // the data

  // NOTE: There's a good chance this whole thing could be sped up a lot by vectorizing this.
  // i.e. change this to a list of row / column vectors
  matrix[numGenes, numCellTypes] sigMat; // the matrix of signature (1000, 2)
}

transformed data
{
  vector<lower=0>[numCellTypes] alpha;

  // not sure this has to be in a loop
  // why not:
  // alpha = 1
  // simply

  // while on the subject - what syntax is this ?

  for (i in 1:numCellTypes)
  {
    alpha[i] = 1;
    // This is an uniform prior over the dirichlet (beta distribution if numCellTypes = 2).
    // Vector of 1s is uniform in the case of  dirichlet distrubution.
  }
}


// The parameters accepted by the model.
parameters
{
  real<lower=0, upper=1> sum_props;
  // Syntax is "simplex[dimensions_of_simplexes] vectorOfSimplexesName[length_of_vector]".
  simplex[numCellTypes] estimatedProportionsVecSimp_unscaled;
  // Syntax is "simplex[dimensions_of_simplexes] vectorOfSimplexesName[length_of_vector]".
  // vector[numCellTypes] estimatedProportionsVecSimp;

  // I should set this up so its a heirarchical prior, and nu is then estimated from the data (see derivation.pdf)
  real<lower=1> nu;

  // Variance parameters, estimated below.
  real<lower=0> sigma;
  real beta0;
}

transformed parameters{
  vector[numCellTypes] estimatedProportionsVecSimp = sum_props*estimatedProportionsVecSimp_unscaled;
}


// The model to be estimated.
model
{
  // alpha is a vector 1's, meaning this dirichlet represents a uniform prior.
  estimatedProportionsVecSimp ~ dirichlet(alpha);

  // Note that the gamma distribution is in terms of shape and rate
  // i.e. y ~ gamma(alpha, beta)
  // this means that for our case: E(nu) = 2/0.01 = 200 !!
  nu ~ gamma(2, 0.1);

  // NOTE: Is it possible to vectorize some of this code????

  // For each genes (1-1000)
  for(geneNum in 1:numGenes)
  {
    real mu;
    mu = 0;

    for(cellTypeNum in 1:numCellTypes)
    {
      mu = mu + sigMat[geneNum, cellTypeNum] * estimatedProportionsVecSimp[cellTypeNum];
    }

    // NOTE: I'm guessing this could be alternatively be implemented as a multivariate liklihood (multivariate t-dist?)?
    // i.e. modelling the correlation across samples ... ?
    // Maybe better if certain assumptions are true, could also hurt if assumptions not met.
    exprMixVec[geneNum] ~ student_t(nu, beta0 + mu, sigma);

  }
}
