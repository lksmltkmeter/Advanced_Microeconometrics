classdef desc
    % desc.m: Module with a few helpful functions for doing simple
    % descriptives
    
    methods (Static)
        
        function binscatter(y, x, nq)
            
            % 1. drop missings
            I = ~isnan(y);
            
            if isdatetime(x)
                I = I &  ~isnat(x);
            else
                I = I & ~isnan(x);
            end
            y = y(I);
            x = x(I);
            
            % 2. find quantiles
            edges = quantile(datenum(x), nq+1);
            x_q = discretize(datenum(x), edges);
            
            % 3. computes means
            if isdatetime(x)
                xx = NaT(nq,1);
            else
                xx = nan(nq,1);
            end
            yy = nan(nq,1);
            for i=1:nq
                I = x_q == i;
                if isdatetime(x)
                    v = round(mean(datenum(x(I))));
                    xx(i) = datetime(v, 'ConvertFrom', 'datenum');
                else
                    xx(i) = mean(x(I));
                end
                yy(i) = mean(y(I));
            end
            
            plot(xx,yy,'-o');
            
        end
        
        function binscatter_by(y, x, cat_var, nq, DOLEGEND)
            cc = unique(cat_var);
            hold on
            for ic=1:numel(cc)
                if iscell(cc(ic))
                    I = strcmp(cat_var, cc{ic});
                else
                    I = cat_var == cc(ic);
                end
                desc.binscatter(y(I), x(I), nq);
            end
            hold off
            
            if nargin < 5 || DOLEGEND
                legend(cc, 'location', 'best');
            end
            
        end
    end
end
