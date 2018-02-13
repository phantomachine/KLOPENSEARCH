function [cq] = cost(q, k, XI)
    % c.m 
    % Provides Cobb-Douglas production's dual cost function
    % XI must be >= 1
    
cq = (q^XI) *(k^(1-XI));

end