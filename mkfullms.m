function [obj, crestfactor, signal]=mkfullms(phases, amplitudes, opts)
    % Generate a full multisine.
    % Input parameters:
    % - *phases*: the value of the phase for each component in radians.
    % - *amplitudes*: the value of the amplitude for each component.
    % - *opts*: struct that contains additional parameters, see mkfullms_opts.
    % Return values:
    % - *obj*: objective that can be used for crest factor minimization, if the amplitudes are not optimized: max(|signal|) 
    % - *crestfactor*: crest factor calculated using peak2rms
    % - *signal*: the generated multisine signal (1D vertical vector of real floating point numbers).
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
