function yhat = mem_tau_exp1(b,X)

% b1 - amplitude of the distribution
% b2 - mean of the distribution
% b3 - kappa of the distribution

yhat = b(1)*exp(-X/b(2)) + b(3);
%Ih =  b(4)*(1./(1+exp(-X/b(5))));
% zi = find(X<b(6));
% Ih(zi) = 0;
%yhat = yhat + Ih;