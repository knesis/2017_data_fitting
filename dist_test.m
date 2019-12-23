function fit = dist_test(probdata,distr)
% dist_test: fits distribution to data and performs chi-squared test
% INPUTS
%   probdata - column vector of 1D data
%   distr - string of desired distribution

fit = {};
fit.pdf = fitdist(probdata, distr);
[fit.h, fit.p, fit.st] = chi2gof(probdata,'CDF', fit.pdf, 'nbins', 9, 'Emin',1);


end