function amplitudes = mkms_default_amplitudes(amplitudes, opts)
    %If the user does not specify the input *amplitudes* (and hence gives []), then initialize them all with 1. If the user specifies them, return the the specified *amplitudes*.
    if isempty(amplitudes)
        amplitudes = ones(opts.n_carriers,1);
    end
