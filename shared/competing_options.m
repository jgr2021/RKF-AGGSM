function options = competing_options(kind)
    if strcmp(kind,'SGaS')
        options = struct('max_iterations',49,'tolerance',1e-2, ...
            'series_terms',31,'series_tolerance',1e-2,'is_particles',1000, ...
            'prior_strength',500,'per_step',false);
    elseif strcmp(kind,'NBKF')
        options = struct('rho',1-exp(-4),'per_step',false);
    else
        error('Unknown competing method');
    end
end
