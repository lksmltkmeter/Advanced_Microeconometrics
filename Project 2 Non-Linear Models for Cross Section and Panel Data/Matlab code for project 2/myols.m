function [b, se, r, cov, R2, s2] = myols(y,x)
% myols estimates the OLS-estimator
% INPUT
%   y:  Dependent variable (Nx1).
%   x:  Explanatory variables (NxK).
% OUTPUT
%   b:  OLS estimates
%   se: Standard errors
%   r:  Vector of residuals
%   cov: Covariance matrix of the OLS coefficients
%   R2: R-squared
%   s2: Residual variance estimate
b=(x'*x)\x'*y;                     % OLS estimates (note the use of the "\"
                                   % operator instead of calculating the
                                   % inverse, i.e. (x'*x)^-1 * x'*y, which
                                   % is less efficient).
r=y-x*b;                           % OLS residuals 
tss=(y-mean(y))'*(y-mean(y));   
rss=r'*r;
ess=tss-rss;                       
R2=ess./tss;                       % R-squared                      
s2=rss/(size(x,1)-size(b,1));      % Residual variance              
cov=s2*(x'*x)^-1;                  % OLS variance-covariance matrix  
se=sqrt(diag(cov));                % OLS standard errors            
return;
