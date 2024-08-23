function X = sample_from_cov(FIM, T)

% This will effectively create "samples" from the passed in FIM about a
% target T. In general, we can use this function to better express the
% geometric relationship between the FIM, Cov ellipse, and physical
% realizations. 

% INPUTS: FIM == 2x2 array; T = 1x2 target.

% sample n points. d is the dimension
n = 100;
d = 2;

% Cov == FIM^-1
Sigma = inv(FIM);

% Generate sample points + center at zero
X = randn(n, d);
X = bsxfun(@minus, X, mean(X));

% inv(chol(cov(X)) normalizes the covariance of X to [1, 0; 0, 1]
X = X * inv(chol(cov(X)));

% Now we overlay our own covariance!
X = X * chol(Sigma);

% Add the target offset and return X
X = X + T;

end

