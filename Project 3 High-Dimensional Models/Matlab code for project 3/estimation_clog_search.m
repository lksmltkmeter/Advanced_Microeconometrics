classdef estimation_clog_search
    %estimation: A matlab class that is really just a collection of different
    %functions put together in one file. 
    
    properties (Constant)
        stepsize_grad = 1e-06; % square root of machine precision [= sqrt(eps)]
        stepsize_hess = 1e-04; % optimal step size is smaller than for gradient
    end
    
    methods (Static)
        
        function [thetahat,se_thetahat,out] = estimate_m(q,theta0,N,opt_search,strCovType) 
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
            if nargin<5;
                strCovType = 'Sandwich'; % default
            end;
            %keyboard
            ll0 = q(theta0);  
            
            %assert(isfinite(sum(ll0)),'Function undefined at initial point. Try different starting values.'); 
            assert(isfinite(nansum(ll0)),'Function undefined at initial point. Try different starting values.'); 
            N = numel(ll0); 
            
            % 1. set up objective function 
            % this must return a scalar value or fminunc will complain 
            Q = @(theta) nansum(q(theta)); 
            
            % 2. estimate parameters 
            %[thetahat,fval,exitflag,out,grad,hess] = fminunc(Q, theta0, options); 
            [thetahat,fval,exitflag,output] = fminsearch(Q,theta0,opt_search) 
                        
            % 3. asymptotic variance matrix 
            % get gradients in N*K vector form
            s = estimation_clog_search.centered_grad(q,thetahat); % s is N*K 
            B = s'*s/N; % s'*s is K*K
            %A = hess/N; % alternatively, "estimation.hessian(q,thetahat)/N"
            keyboard
            A = estimation_clog_search.hessian(q,thetahat)/N; 
            switch strCovType
                case 'Sandwich'
                    cov = 1/N * (A\B/A);
                case 'Outer Product'
                    cov = 1/N * B^-1;
                case 'Hessian'
                    cov = 1/N * A^-1;
                otherwise 
                    error('Unexpected value of input strCovType, %s. Must be ''Sandwich'', ''Outer Product'', or ''Hessian''.',strCovType); 
            end; 
            
            % 4. standard errors
            se_thetahat = sqrt(diag(cov)); 
            
            % 5. print results 
            estimation_clog_search.print_results(thetahat,se_thetahat); 
            
        end; 
        
        function [thetahat,se] = estimate_model(model,y,x,options,strCovType)
            [N,K] = size(x); 
            theta0 = model.starting_values(y,x); 
            q = @(theta) model.criterion(y,x,theta); 
            [thetahat,se] = estimation_clog_search.estimate_m(q,theta0,N,options,strCovType); 
        end; 
        
        function [thetahat,se_thetahat,out] = estimate_bhhh(model,y,x,options,theta0)
            q = @(theta) model.criterion(y,x,theta); 
            Q = @(theta) estimation_clog_search.crit_for_bhhh(q,theta); 
            
            if nargin<5; 
                theta0 = model.starting_values(y,x); 
            end; 
            
            options.Algorithm = 'trust-region'; 
            options.GradObj = 'on'; 
            options.Hessian = 'on'; 
            [thetahat,fval,exitflag,out,grad,hess] = fminunc(Q, theta0, options); 
            
            ll = q(thetahat); 
            N = numel(ll); 
            cov = hess^-1; 
            se_thetahat = sqrt(diag(cov)); 
            
            estimation_clog_search.print_results(thetahat,se_thetahat); 
        end; 
        
        function [ Q,g,h ] = crit_for_bhhh(q,theta)
            Q = sum(q(theta)); 
            s = estimation_clog_search.centered_grad(q,theta); 
            h = s'*s; 
            g = sum(s,1); 
        end; 

        function [ grad ] = centered_grad( fun , x0 , h )
            % Computes an estimate of the gradients of a vectorized function using
            % centered differences and a relativestep size of 1e-08. If x0 is zero at
            % any element, it will just take a step of 1e-08.
            
            f0 = fun(x0);
            N = numel(f0);
            grad = nan(N,numel(x0));
            
            if nargin<3; 
                h = estimation_clog_search.stepsize_grad; 
            end; 
            
            for i=1:numel(x0);
                
                % 1. find forward and backward step 
                x1 = x0;  % forward step
                x_1 = x0; % backwards step
                if x0(i)~=0;
                    x1 = x0;
                    x1(i) = x0(i)*(1+h);
                    x_1(i) = x0(i)*(1-h);
                else % if x0(i)==0, then we cannot take a relative step, so we take an absolute step instead
                    x1(i) = h;
                    x_1(i) = -h;
                end;
                steplength = x1(i)-x_1(i);
                
                % 2. compute newton quotient 
                grad(:,i) = (fun(x1)-fun(x_1))/steplength;
            end;
            
        end
        
        function hess = hessian( fhandle , x0 )
            % Computes the hessian of the input function at the point x0 
            
            % Initialization
            K = size(x0,1);
            f2 = zeros(K,K); % double step
            f1 = zeros(K,1); % single step
            h_rel = estimation_clog_search.stepsize_hess; 
                        
            % Step size 
            dh = h_rel*x0;      % use relative steps
            dh(x0==0) = h_rel;  % ... except when x0(k)==0
            
            % Initial point 
            f0vec = fhandle(x0); 
            f0 = sum(f0vec,1);
            
            % Single forward step 
            for k=1:K; 
                x1 = x0;% + ee(:,k);
                x1(k) = x0(k) + dh(k); 
                f1(k,1) = sum(fhandle(x1),1);
            end;
            
            % Double forward step 
            for k=1:K; 
                for j=1:k; % only loop to the diagonal
                    
                    % 1. find the new point (after double-stepping) 
                    x2 = x0; 
                    if k==j; % for numerical stability we do this in a single step
                        x2(k) = x0(k) + dh(k) + dh(k); 
                    else % 
                        x2(k) = x0(k) + dh(k); 
                        x2(j) = x0(j) + dh(j); 
                    end; 

                    % 2. compute function value 
                    f2(k,j) = sum(fhandle(x2),1);
                    
                    % 3. fill out above the diagonal 
                    if k ~= j; 
                        f2(j,k) = f2(k,j);
                    end;
                end;
            end;
            
            f1rep = repmat(f1,[1,K]); 
            
            % F1(k,j) gives the effect of first taking a step in direction
            % k and then a step in direction j, purely from the first
            % derivative (f1). 
            F1 = f1rep + f1rep'; 
            
            % denominator is K*K outer product of the steps
            hess = (f2 - (F1 - f0)) ./ (dh*dh');
        end; 
        
        function print_results(theta,se)
            for i=1:numel(theta); 
                if nargin==1; 
                    fprintf('theta_%i: %8.4f\n',i,theta(i)); 
                else % also std.errs. 
                    fprintf('theta_%i: %8.4f  (se = %6.4f, t = %5.2f)\n',i,theta(i),se(i),theta(i)/se(i)); 
                end; 
            end; 
        end;         
        
        function [xx,yy] = plot_data_and_prediction_vs_x2(y,x,model,theta)
            % Shows x(:,2) against y and predicted y
            
            % Compute fitted values from the model 
            yhat = model.predict(x,theta); 
            
            % Plot it 
            A = sortrows([x(:,2),yhat]); 
            xx = A(:,1); 
            yy = A(:,2); 
            p = plot(x(:,2),y,'.',xx,yy,'-'); 
            legend('Data','Prediction','Location','Best');
            set(gcf,'color','w');
            set(gca,'box','off');
            set(gca,'fontsize',14);
            p(1).MarkerSize = 20; 
            p(2).LineWidth = 3; 
        end; 
        
        function rmse = compute_rmse(y,x,model,theta)
            yhat = model.predict(x,theta);
            rmse = sqrt(mean((y-yhat).^2));
            fprintf('RMSE = %f.\n',rmse);
        end;
        
        function [xx,ff] = plot_numerical_noise(q,theta,iPar,h)
            % This function plots the criterion function, sum(q), around
            % the given parameter vector, theta, where only theta(iPar) is
            % changed. The function is evaluated 100 times between the
            % value theta(iPar)*(1-h) to theta(iPar)*(1+h). By setting h=1e-08, for
            % example, one can explore how much numerical noise is in the
            % criterion function. 
            
            nEval = 100; 
            
            theta2 = theta(:); % useful if theta is a matrix

            collapse = @(vec) sum(vec);
            
            % Set up the function handle so that all parameters except
            % number iPar are kept at the value in theta, while theta(iPar)
            % is changed. 
            if iPar==1; % only first parameter is changed
               fh = @(t1) collapse(q([t1;theta2(2:end)]));
            elseif iPar<numel(theta2); % change a parameter in the middle
                fh = @(t1) collapse(q([theta2(1:iPar-1);t1;theta2(iPar+1:end)]));
            else % only last parameter is changed 
               fh = @(t1) collapse(q([theta2(1:end-1);t1]));
            end; 
            
            % Values of theta
            xx = linspace(theta2(iPar)*(1-h),theta2(iPar)*(1+h),nEval)';
            % Function values 
            ff = nan(nEval,1);
            for p=1:nEval;
                ff(p) = fh(xx(p)); % evaluate the function 
            end;
            
            plot(xx,ff);
            
            ylabel('$Q(\theta)$','interpreter','latex'); 
            xlabel(sprintf('$\\theta_{%i}$',iPar),'interpreter','latex'); 
            set(gcf,'color','w');set(gca,'box','off'); set(gca,'fontsize',14);
        end;

    end
    
end

