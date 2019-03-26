function [mu, sigma2] = posteriorGaussian(mu_x, sigma2_x, y, sigma2_r)
%posteriorGaussian performs a single scalar measurement update with a
%measurement model which is simply "y = x + noise".
%
%Input
%   MU_P            The mean of the (Gaussian) prior density.
%   SIGMA2_P        The variance of the (Gaussian) prior density.
%   SIGMA2_R        The variance of the measurement noise.
%   Y               The given measurement.
%
%Output
%   MU              The mean of the (Gaussian) posterior distribution
%   SIGMA2          The variance of the (Gaussian) posterior distribution

%Your code here

end

% A few random test vectors
n = 10;
mu_p = 5*randn(1,n);
sigma2_p = 3*rand(1,n);
y = mu_p+5*randn(1,n);
sigma2_r = 3*rand(1,n);
tol = 1e-8;

for k = 1:n
    [mu1,sigma1] = posteriorGaussian(mu_p(k), sigma2_p(k), y(k), sigma2_r(k));
    [mu2,sigma2] = reference.posteriorGaussian(mu_p(k), sigma2_p(k), y(k), sigma2_r(k));
    assert(abs(mu1-mu2) < tol, 'mean is not correct');
    assert(abs(sigma1-sigma2) < tol, 'variance is not correct');
end


%% Fixed test case with public true values
mu_p = 3.3575;
sigma2_p = 2.2732;
sigma2_r = 3.4704;
y = 5.4159;
mu2 = 4.1722;
sigma2 = 1.3735;

[mu1,sigma1] = posteriorGaussian(mu_p, sigma2_p, y, sigma2_r);
tol = 1e-4;

assert(abs(mu1-mu2) < tol, 'mean is not correct');
assert(abs(sigma1-sigma2) < tol, 'variance is not correct');

