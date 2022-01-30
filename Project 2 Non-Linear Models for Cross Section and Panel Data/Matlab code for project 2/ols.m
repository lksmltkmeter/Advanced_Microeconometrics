classdef ols
    %OLS functions for estimation, data generation, etc.

    properties
    end

    methods (Static)

        function [y,x] = simulate_dataset(N,theta)
            K = numel(theta)-1;

            % 0. unpack parameters
            beta = theta(1:K);
            sigma = theta(K+1);

            % 1. simulate x variables.
            % x has a constant and one variable that distributed standard
            % normal.
            x = [ones(N,1),randn(N,K-1)];

            % 2. simulate y values
            % Add a normal error term with std.dev. = sigma.
            y = x*beta + sigma * randn(N,1);
        end

        function theta0 = starting_values(y,x)
            % Does not really make sense here but we have to fill it out
            theta0 = ols.estimate(y,x);
        end

        function q = criterion(y,x)
            q = nls.criterion(y,x);
        end

        function [thetahat,theta_se] = estimate(y,x,solverOptions)
            % Note, the input "solverOptions" is not actually needed for
            % OLS, since we find the estimator in closed form.
            N = size(x,1);
            K = size(x,2);
            % check inputs
            if size(y,1)~=N
                error('Unexpected dimension of y. Should be %i*1, got %i*%i.',N,size(y,1),size(y,2));
            end
            
            % 1. estimate beta
            % Note, for matrices A and B: A\B is the same as A^-1 * B, but it is faster for the computer
            betahat = (x'*x)\x'*y;

            % 2. estimate sigma
            res = y - x*betahat;
            sigma2hat = 1/(N-K) * (res'*res);

            % 3. standard errors
            cov = sigma2hat * inv(x'*x);
            se = sqrt(diag(cov));
            se_beta = se(1:K);

            thetahat = [betahat;sqrt(sigma2hat)];
            theta_se = [se_beta;nan]; % we have not filled out the std.err. for the variance parameter
        end

    end

end

