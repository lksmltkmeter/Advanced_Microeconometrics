classdef tobit_1_logistic
    %tobit
    %pd = makedist('Logistic')
    %log_cdf = cdf(pd)
    %log_cdf(pd,xb/sigma)

    properties
        % No properties
    end

    methods (Static)

        function [y,x] = simulate_dataset(N,theta)
            % 0. unpack parameters from theta
            beta = theta(1:end-1);
            sigma = theta(end);

            % 1. simulate x variables.
            x = [ones(N,1),randn(N,1)];

            % 2. simulate y values
            ystar = x*beta + sigma*randn(N,1);
            y = max(ystar,0);
        end

        function theta0 = starting_values(y,x)
            idx = y>0;
            theta0 = ols.estimate(y(idx),x(idx,:));
            theta0(end) = theta0(end)*1.5;  % recommended to start with a larger sigma
        end

        function yhat = predict(x,theta)
            beta = theta(1:end-1);
            sigma = theta(end);

            xb = x*beta;

            yhat = xb.*normcdf(xb/sigma) + sigma * normpdf(xb/sigma);
        end

        function q = criterion(y,x,theta)
            % 0. unpack
            beta = theta(1:end-1);
            sigma = theta(end);
            pd = makedist('Logistic');
          
            % 2. useful intermediate variables
            xb = x*beta;
            res = y - xb;

            % 3. criterion
            ll = (y==0).* log(1 - cdf(pd,xb/sigma)) ...
               + (y>0) .* (log((1/sigma)*pdf(pd,res/sigma)));
            q = -ll;
        end


    end

end

