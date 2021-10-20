function [chi22]=chitestquantities(paramt,parame,Ce)

% chitestquantities calculates the chi2 increment for a 2 by 1 vector and
% its covariance matrix
% paramt    true parameters (from simulation)
% parame    estimated parameters (from reconstruction)
% Ce        2 by 2 covariance matrix

delta=parame-paramt;
chi22=delta'*inv(Ce)*delta;