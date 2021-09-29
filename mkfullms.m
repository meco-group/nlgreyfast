function [obj, crestfactor, signal]=mkfullms(phases, amplitudes, opts)
    % generate a full multisine
    amplitudes = mkms_default_amplitudes(amplitudes, opts);
    assert(all(size(phases)==[opts.n_carriers 1]))
    assert(all(size(amplitudes)==[opts.n_carriers 1]))
    if opts.type == mkms_type.full % ~isfield(opts,'type') | opts.type == [] %this is an untested change
        complex_values=[zeros(opts.n_leading_zeros,1); amplitudes.*exp(1j*phases); zeros(opts.n_trailing_zeros,1)];
    elseif opts.type == mkms_type.custom
        complex_values=zeros(opts.n_samples,1);
        complex_values(opts.fft_interesting_indexes) = amplitudes.*exp(1j*phases);
    end
    signal=2*real(ifft(complex_values));
    obj=max(abs(signal));
    if nargout>=2
        crestfactor=peak2rms(signal);
    end
end
