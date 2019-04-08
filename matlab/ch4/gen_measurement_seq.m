function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
% sequence X using a linear measurement model. Measurement noise is assumed to be 
% zero mean and Gaussian.
%
%Input:
%   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
%   H           [m x n] Measurement matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% your code here
end

absTol = 1e-1;
relTol = 5e-2;

N = 50000;

n = randi(1,1);
m = randi(n,1);

% Define state sequence
X = zeros(n,N+1);

% Define measurement model
H = 1;
R = .5^2;

% Generate measurements
Y = genLinearMeasurementSequence(X, H, R);

% PLot results
figure(1);clf;hold on;
plot(0:10,X(1,1:11), '--k');
plot(1:10, Y(1,1:10), '*r');
legend('State sequence', 'Measurements')
title('Your solution');
xlabel('k');
ylabel('position');

assert(size(Y,1) == m, 'Y has the wrong measurement dimension');
assert(size(Y,2) == N, 'Y should have N columns');

Rest = cov((Y-H*X(:, 2:N+1))');
assert(all(all((mean(Y-H*X(:, 2:N+1)) < absTol))), 'Measurement noise is not zeros mean');
assert(all(all((abs(Rest-R) < relTol*R))), 'Measurement noise covariance is not within tolerances');

%%%%%%%%%%%%

absTol = 1e-1;
relTol = 5e-2;


N = 10000;

n = randi(5,1);
m = randi(n,1);

N = m*N;

% Define state sequence
X = rand(n,N+1);

% Define measurement model
H = randn(m,n);
V = rand(m,m);
R = V*diag(10*rand(m,1))*V';

% Generate measurements
Y = genLinearMeasurementSequence(X, H, R);

Rest = cov((Y-H*X(:, 2:N+1))');
assert(size(Y,1) == m, 'Y has the wrong measurement dimension');
assert(size(Y,2) == N, 'Y should have N columns');

Rest = cov((Y-H*X(:, 2:N+1))');
assert(all(all((mean(Y-H*X(:, 2:N+1),2) < absTol))), 'Measurement noise is not zeros mean');
assert(all(all((abs(Rest-R) < relTol*R))), 'Measurement noise covariance is not within tolerances');


