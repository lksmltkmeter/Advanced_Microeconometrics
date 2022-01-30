classdef clad
    %clad

    properties
        % No properties
    end

    methods (Static)


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

            % 2. useful intermediate variables
            xb = x*beta;
            res = y - xb;

            % 3. criterion
            q = abs(y-max(xb,0));
        end


    end

end

