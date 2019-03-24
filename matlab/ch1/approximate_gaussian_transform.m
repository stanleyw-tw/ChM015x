function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f


if nargin < 4
    N = 5000;
end

%Your code here

end

run visible_shared_variables
% transform Gaussian distribution and calculate mean and covariance by sampling
[mu_y, Sigma_y, xs] = approxGaussianTransform(mu_x, Sigma_x, f);

[mu_ref, Sigma_ref] = reference.approxGaussianTransform(mu_x, Sigma_x, f);
dm = max(abs(mu_ref - mu_y));
assert(all(abs(mu_ref - mu_y) < 0.1), 'mean is outside of tolerance by %f', dm)
ds = max(abs(Sigma_ref(:) - Sigma_y(:)));
assert(all(abs(Sigma_ref(:) - Sigma_y(:)) < 0.2), 'covariance is outside of tolerance by %f', ds)


run visible_shared_variables
% transform Gaussian distribution and calculate mean and covariance by sampling
N = 1;
[mu_y, Sigma_y, ys] = approxGaussianTransform(mu_x, Sigma_x, f, N);

dElems = diag(Sigma_y);
assert(all(~(dElems < 1e-10)), 'Covariance should not be zero when using only one sample point')


run visible_shared_variables
% transform Gaussian distribution and calculate mean and covariance by sampling
N = 10+round(rand*1000);
[mu_y, Sigma_y, ys] = approxGaussianTransform(mu_x, Sigma_x, f, N);

assert(size(ys,2) == N, 'Number of samples used and returned should be controlled by the fourth input parameter')


