classdef poisson
    %tobit

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

            yhat_Poisson = exp(x*beta);

        end



        function q = criterion(y,x,theta)
            % 0. unpack
            beta = theta(1:end-1);
            sigma = theta(end);

            % 2. useful intermediate variables
            xb = x*beta;
            %res = y - xb;

            % 3. criterion
            %ll = (y==0).* log(1 - normcdf(xb/sigma)) ...
               %+ (y>0) .* (-.5*log(2*pi) - log(sigma) - .5*res.^2/(sigma^2) );
            ll= -exp(xb)+y.*xb;
            q = -(ll);
        end


    end

end

