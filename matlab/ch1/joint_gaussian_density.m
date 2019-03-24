function [mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)
%jointGaussian calculates the joint Gaussian density as defined
%in problem 1.3a. 
%
%Input
%   MU_X        Expected value of x
%   SIGMA2_X    Covariance of x
%   SIGMA2_R    Covariance of the noise r
%
%Output
%   MU          Mean of joint density 
%   SIGMA       Covariance of joint density


%Your code here

end


%%%parameters from problem
% parameters as given in problem
mu_x = 19;
sigma2_x = 5^2;
sigma2_r = 2^2;

% tolerance
tol = 1e-8;

[mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r);
[mu_ref, Sigma_ref] = reference.jointGaussian(mu_x, sigma2_x, sigma2_r);
assert(all(abs(mu-mu_ref) < tol), 'mean is wrong');
assert(all(abs(Sigma(:)-Sigma_ref(:)) < tol), 'covariance is wrong');

%%%random test parameters
% random test case 
mu_x1 = 5+randn*10;
P_x1 = (randn*3)^2;
sigma2_r1 = (randn*2)^2;

% tolerance
tol = 1e-8;

[mu, Sigma] = jointGaussian(mu_x1, P_x1, sigma2_r1);
[mu_ref, Sigma_ref] = reference.jointGaussian(mu_x1, P_x1, sigma2_r1);
assert(all(abs(mu-mu_ref) < tol), 'mean is wrong');
assert(all(abs(Sigma(:)-Sigma_ref(:)) < tol), 'covariance is wrong');

%%% Check the dimensionality of the output 
mu_x = 19;
sigma2_x = 5^2;
sigma2_r = 2^2;

[mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r);
assert(size(mu,1) == 2, 'mu should have two rows, when all inputs are scalar');
assert(size(mu,2) == 1, 'mu should have one column, when all inputs are scalar');
assert(size(Sigma,1) == 2, 'Sigma should have two rows, when all inputs are scalar');
assert(size(Sigma,2) == 2, 'Sigma should have two column, when all inputs are scalar');

