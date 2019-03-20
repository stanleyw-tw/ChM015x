function [mu_y, Sigma_y] = affineGaussianTransform(mu_x, Sigma_x, A, b)
%affineTransformGauss calculates the mean and covariance of y, the 
%transformed variable, exactly when the function, f, is defined as 
%y = f(x) = Ax + b, where A is a matrix, b is a vector of the same 
%dimensions as y, and x is a Gaussian random variable.
%
%Input
%   MU_X        [n x 1] Expected value of x.
%   SIGMA_X     [n x n] Covariance of x.
%   A           [m x n] Linear transform matrix.
%   B           [m x 1] Constant part of the affine transformation.
%
%Output
%   MU_Y        [m x 1] Expected value of y.
%   SIGMA_Y     [m x m] Covariance of y.

%Your code here
    mu_y = A*mu_x+b;
    Sigma_y = A*Sigma_x*A';
end


%% compare resulting mean and cov to reference 2-D
A3 = ...
  [4 1;
   -1 1];

P3 = ...
  [4 2;
   2 5];

b3 = [-0; 3];

mu3 = [3; -6];

s3 = ...
  [0 2;
   2 1];

x_cov3 = ...
  [85 -5;
   -5 5];

x_mean3 = [6; -6];

tol = 1e-5;

[x_mean, x_cov] = affineGaussianTransform(mu3, P3, A3, b3);
assert(all(abs(x_mean3 - x_mean) < tol), 'mean is outside of tolerance by %f', max(abs(x_mean3 - x_mean)))
assert(all(abs(x_cov3(:) - x_cov(:)) < tol), 'covariance is outside of tolerance by %f', max(abs(x_cov3(:) - x_cov(:))))



run('visible_shared_variables')
% transform Gaussian distribution and calculate mean and covariance
[x_mean, x_cov] = affineGaussianTransform(mu1, P1, A1, b1);
[ref_mean, ref_cov] = reference.affineGaussianTransform(mu1, P1, A1, b1);
assert(all(abs(ref_mean - x_mean) < tol), 'mean is outside of tolerance by %f', max(abs(ref_mean - x_mean)))
assert(all(abs(ref_cov(:) - x_cov(:)) < tol), 'covariance is outside of tolerance by %f', max(abs(ref_cov(:) - x_cov(:))))


run('visible_shared_variables')
% transform Gaussian distribution and calculate mean and covariance
[x_mean, x_cov] = affineGaussianTransform(mu2, P2, A2, b2);
[ref_mean, ref_cov] = reference.affineGaussianTransform(mu2, P2, A2, b2);
assert(all(abs(ref_mean - x_mean) < tol), 'mean is outside of tolerance by %f', max(abs(ref_mean - x_mean)))
assert(all(abs(ref_cov(:) - x_cov(:)) < tol), 'covariance is outside of tolerance by %f', max(abs(ref_cov(:) - x_cov(:))))


run('visible_shared_variables')
% transform Gaussian distribution and calculate mean and covariance
[x_mean, x_cov] = affineGaussianTransform(mu5, P5, A5, b5);
[ref_mean, ref_cov] = reference.affineGaussianTransform(mu5, P5, A5, b5);
assert(all(abs(ref_mean - x_mean) < tol), 'mean is outside of tolerance by %f', max(abs(ref_mean - x_mean)))
assert(all(abs(ref_cov(:) - x_cov(:)) < tol), 'covariance is outside of tolerance by %f', max(abs(ref_cov(:) - x_cov(:))))


