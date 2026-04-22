function x = copnorm(x)
% COPNORM Copula normalisation
%   cx = copnorm(x) returns standard normal samples with the same empirical
%   CDF value as the input. Operates along the first axis.
%   Equivalent to cx = norminv(ctransform(x))
%
[~,x] = sort(x, 1);
[~,x] = sort(x, 1);
x = x / (size(x, 1) + 1);
x = -sqrt(2).*erfcinv(2*x);

%% check if there are too many  repeated ranks
D = size(x,2); %variables' amount
for d = 1:D
    unique_vals = numel(unique(x(:,d)));
    redundancy = 1 - unique_vals / size(x, 1);
    if redundancy > 0.1
        warning('[copnorm] Column %d has %.2f%% repeated copula ranks.', d, redundancy * 100);
    end
end