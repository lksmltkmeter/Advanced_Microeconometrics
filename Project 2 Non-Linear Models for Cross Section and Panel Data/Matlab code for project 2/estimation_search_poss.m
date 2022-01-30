classdef estimation_search_poss
    %estimation: A matlab class that is really just a collection of different
    %functions put together in one file.

    properties
    end

    methods (Static)

        function [thetahat,se_thetahat, cov] = estimate_m(q,theta0,N,opt_search,strCovType)
           
            % for fminunc use estimate_m(q,theta0,N,opt,strCovType)
            % Estimate parameters of an M-estimator
            %
            % INPUTS:
            %   q:          The criterion function (must return a vector)
            %   theta0:     Starting values for the optimizer
            %   options:    Options for the gradient-based optimizer, fminunc
            %   strCovType: [optional] If provided, must be 'Outer Product'
            %                or 'Sandwich' (default).
            % OUTPUTS:
            %   thetahat:   Parameter estimate
            %   se_thetahat: Standard errors
            if nargin<5
                strCovType = 'Sandwich'; % default
            end

            % 1. set up objective function
            % this must return a scalar value or fminunc will complain
            Q = @(theta) sum(q(theta));

            % 2. estimate parameters
            %[thetahat,fval,exitflag,out,grad,hess] = fminunc(Q, theta0, opt);
            %[thetahat,fval,exitflag,output, grad,hess] = fminsearch(Q,theta0,opt_search)
            [thetahat,fval,exitflag,output] = fminsearch(Q,theta0,opt_search)


            % 3. asymptotic variance matrix
            % get gradients in N*K vector form
            %thetahat_poss = thetahat(1:10,:); 
            s = estimation_search.centered_grad(q,thetahat); % s is N*K
            B = s'*s/N; % s'*s is K*K
            %for fminunc
            A = hessian(Q,thetahat)/N;
            
            %A = hess/N; % alternatively, "estimation.hessian(q,thetahat)/N"
            switch strCovType
                case 'Sandwich'
                    cov = 1/N * (A\B/A);
                case 'Outer Product'
                    cov = 1/N * B(1:13,1:13)^-1;
                 case 'Hessian'
                    cov = 1/N * A^-1;
                otherwise
                    error('Unexpected value of input strCovType, %s. Must be ''Sandwich'', ''Outer Product'', or ''Hessian''.',strCovType);
            end
            
            % 4. standard errors
            se_thetahat = sqrt(diag(cov));

            % 5. print results
            %estimation_search_poss.print_results(thetahat,se_thetahat);

        end

        function grad = forward_diff( fun,x0,h )
            f0 = fun(x0);
            N = numel(f0);
            grad = nan(N,numel(x0));
            if nargin<3
                h = 1.49e-08; % square root of machine precision [= sqrt(eps)]
            end

            f0 = fun(x0);
            for i=1:numel(x0)
                x1 = x0;  % forward step
                if x0(i)~=0
                    x1 = x0;
                    x1(i) = x0(i)*(1+h);
                else % if x0(i)==0, then we cannot take a relative step, so we take an absolute step instead
                    x1(i) = h;
                end
                step = x1(i) - x0(i);
                grad(:,i) = (fun(x1)-f0)/step;
            end

        end

        function [ grad ] = centered_grad( fun , x0 )
            % Computes an estimate of the gradients of a vectorized function using
            % centered differences and a relativestep size of 1e-08. If x0 is zero at
            % any element, it will just take a step of 1e-08.

            f0 = fun(x0);
            N = numel(f0);
            grad = nan(N,numel(x0));
            h = 1.49e-08; % square root of machine precision [= sqrt(eps)]

            for i=1:numel(x0)
                x1 = x0;  % forward step
                x_1 = x0; % backwards step
                if x0(i)~=0
                    x1 = x0;
                    x1(i) = x0(i)*(1+h);
                    x_1(i) = x0(i)*(1-h);
                else % if x0(i)==0, then we cannot take a relative step, so we take an absolute step instead
                    x1(i) = h;
                    x_1(i) = -h;
                end
                step = x1(i)-x_1(i);
                grad(:,i) = (fun(x1)-fun(x_1))/step;
            end

        end

        function print_results(theta,se)
            for i=1:numel(theta-1)
                if nargin==1
                    fprintf('theta_%i: %8.4f\n',i,theta(i));
                else % also std.errs.
                    fprintf('theta_%i: %8.4f  (se = %6.4f, t = %5.2f)\n',i,theta(i),se(i),theta(i)/se(i));
                end
            end
        end
    end

end

