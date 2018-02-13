function out = cpi(Ph,Pf,EPSILON,THETA)

    % CES aggregator
    
    % THETA - Home bias share
    % EPSILON - s.t. 1/EPSILON - 1 is elasticity of substition btw x-y
    
    if EPSILON == Inf
        out = (Ph^THETA)*(Pf^(1-THETA))/((THETA^THETA)*((1-THETA)^(1-THETA)));
    elseif EPSILON == 1
        warning('cpi:Error','Limit EPSILON = 1 not admissible!')
    elseif EPSILON == 0
        warning('cpi:Error','Leotief EPSILON = 0 not admissible!')
    elseif EPSILON ~= Inf || EPSILON ~= 0 || EPSILON ~= 1
        out = ( (THETA^(EPSILON/(EPSILON-1)))*(Ph)^(1/(1-EPSILON)) ...
                + ((1-THETA)^(EPSILON/(EPSILON-1)))...
                                *(Pf)^(1/(1-EPSILON)) )^(1-EPSILON);
    end