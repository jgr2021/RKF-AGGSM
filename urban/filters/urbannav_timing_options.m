function options = urbannav_timing_options(options)
% Optional measurement controls; legacy calls retain their loss output.
% Timers exclude initialization and calculate_loss. Per-step timing is a
% separate diagnostic pass; sequence benchmarks leave it disabled.
    defaults = struct('enabled', false, 'per_step', false, 'compute_loss', true);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(options, names{k})
            options.(names{k}) = defaults.(names{k});
        end
        validateattributes(options.(names{k}), {'logical'}, {'scalar'});
    end
    options.enabled = options.enabled || options.per_step;
end
