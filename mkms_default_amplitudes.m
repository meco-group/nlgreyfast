function amplitudes = mkms_default_amplitudes(amplitudes, opts)
    %if input amplitudes is [] return ones, otherwise return the input amplitudes
    if isempty(amplitudes)
        amplitudes = ones(opts.n_carriers,1);
    end
