function opts = optdefaults(opts, defaults)
    % opts = optdefaults(opts, defaults)
    %
    %   Set default field values in opts array.  For any field of defaults
    %   not set in opts, adds that field/value to opts.
    %
    
    names = fieldnames(defaults);
    
    for i=1:length(names)
        name = names{i};
        
        if ~isfield(opts, name)
            opts.(name) = [defaults.(name)];
        end
    end
end