function print_par(theta, se, lab, header)
% INPUTS
%   theta: K-vector of parameters
%   se: K-vector of standard errors 
%   lab: cell array of K names of variables 
%   header: string, the title of the table 
% [no outputs]
    K = numel(theta); 
    assert(numel(se) == K, 'theta and se do not have the same number of elements'); 
    assert(numel(lab) == K, 'theta and lab do not have the same number of elemetns'); 
    
    fprintf('--- %s ---\n', header); 
    for k=1:numel(lab)
        fprintf('%20s: %12.4g (t=%8.2f)\n', lab{k}, theta(k), theta(k)/se(k)); 
    end
end
