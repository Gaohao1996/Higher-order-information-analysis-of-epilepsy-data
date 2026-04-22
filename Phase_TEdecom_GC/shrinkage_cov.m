function Sigma_shrink = shrinkage_cov(X)

    [N, d] = size(X);
    mu = mean(X, 1);
    Xc = X - mu;
    S = (Xc' * Xc) / (N - 1);  % 

    % Shrinkage target: scaled identity matrix
    var_avg = trace(S) / d;
    T = var_avg * eye(d);

    % Ledoit-Wolf shrinkage coefficient
    diff = Xc.^2 - var(Xc, 1);
    beta_hat = sum(diff(:).^2) / N;
    alpha_hat = norm(S - T, 'fro')^2;
    lambda = min(max(beta_hat / alpha_hat, 0), 1);  % Clamp to [0,1]

    % Shrinkage estimator
    Sigma_shrink = (1 - lambda) * S + lambda * T;
end
