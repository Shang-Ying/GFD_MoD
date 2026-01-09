% Unconditional exponential covariogram model

function [cov_uu]=fct_UnCovMat(pos_est,sill,range,nugget)
Duu = pdist2(pos_est,pos_est); % unknown to unknown distance matrix
cov_uu = nugget+sill*exp(-1*Duu/range); % unknown to unknown covariance matrix
