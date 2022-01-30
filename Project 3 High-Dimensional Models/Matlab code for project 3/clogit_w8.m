classdef clogit_w8
    %logit
    
    properties (Constant) 
        % not used 
    end
    
    methods (Static) 
        
        function [y,x] = simulate_dataset(N,J,theta)
            % 0. unpack parameters from theta
            beta = theta; 
            K = numel(beta); 
            
            % 1. simulate x variables (characteristics).
            x = randn(N,J,K); %simulate x as NxJxK random normal
            
            % 2. utility
            x_r = reshape(x,[N*J, K]);  % JN*1 to prepare for product 
            V = reshape(x_r*theta,N,J); % N*J
            E = gevinv(rand(N,J));  
             
            U = V + E ; % N*J
            
            % 3. chosen alternative maximizes utility 
            [maxU,y] = max(U,[],2);       
            

        end
        
        function [ccp,logccp] = choice_probabilities(x,theta)
            [N,J,K] = size(x); 
            assert(K~=1,'Expecting x to be a 3-dimensional array.'); 
            x_r = reshape(x,N*J,K); 
            %keyboard
            % 1. utilities (J*1 vector) 
            
            u = x_r*theta; 
            utils = reshape(u,N,J); 
            
            % 2. subtract maximum for numerical stability 
            maxU = max(utils,[],2);       % N*1
            maxU_N_J = maxU(:,ones(1,J)); % N*J
            utils = utils - maxU_N_J; 
            
            % 3. denominator 
            denom = nansum(exp(utils),2);      % N*1
            %denom = sum(exp(utils),2);      % N*1
            denom_N_J = denom(:,ones(1,J)); % N*J 
            %keyboard
            % 4. conditional choice probabilities (CCPs), J*1 vector
            ccp = exp(utils) ./ denom_N_J;  % N*J
            logccp = utils - log(denom_N_J);% N*J 
        end 

        function theta0 = starting_values(y,x)
            if numel(size(x))>=3
                [N,J,K] = size(x); 
                theta0 = 0.01 * ones(K, 1); 
            else
                [N,K] = size(x); 
                theta0 = 0.01 * ones(K, 1); 
            end
            %theta0 = zeros(K,1); 
        end 

        function q = criterion(y,x,theta)
             % Returns the _negative_ log likelihood
            % INPUTS: 
            %   y,x:    The data
            %   theta:  Parameters
            % OUTPUTS
            %   q:      The criterion function (in vector form)
            
            [N,J,K] = size(x); 
            
            %y_ = reshape(y,N*J,1);
            [ccp,logccp] = clogit_w8.choice_probabilities(x,theta); 
            %keyboard
            %idx_choice = sub2ind(size(ccp),(1:N)',y);
            %reshape ccp 
            %y_ = reshape(y, N,J);
            %idx_choice = sub2ind(size(ccp),(1:N)',y_); 
            ccp_r = reshape(ccp, N*J,1);
            %keyboard
            %SET ALL NAN TO 0
            %ccp_r(isnan(ccp_r))=0; 
            %y(isnan(y))=0; 
            %idx_choice = sub2ind(size(ccp),(1:N*J)',y_);
            %q = -logccp(idx_choice); 
           
            q = (y - ccp_r).^2;
            %q = reshape(q,N,J); % recall that N = N*T so q should be reshaped to N*TxJ
            %q = - (y - ccp_r(idx_choice)).^2;
            %q=nansum(q,2); %sum over rows (see dicussion forum suggestion by Johan & Christian)
        end        
        
    end
    
end

