function opts = mkfullms_opts(n_samples,fs,fmin,fmax)
    % Create options structure for mkfullms, in order to 
    % generate a full multisine with a total number of samples *n_samples*, sinusoidal component frequencies between *fmin* and *fmax*, sampled at *fs*.
    opts.df = fs/n_samples; %frequency difference between two consecutive carrier frequencies
    opts.n_leading_zeros = floor(fmin/opts.df);
    opts.n_carriers = round((fmax-fmin+0.5*opts.df)/opts.df);
    opts.n_trailing_zeros = n_samples-opts.n_leading_zeros-opts.n_carriers; 
    opts.n_samples = n_samples; %number of samples in time domain signal
    opts.fs = fs; %sampling frequency
    opts.f = ((opts.n_leading_zeros-1+(1:(opts.n_carriers)))*opts.df).'; %all carrier frequencies
    opts.fmin = min(opts.f); %minimum carrier frequency
    opts.fmax = max(opts.f); %maximum carrier frequency
    opts.fft_interesting_indexes=opts.n_leading_zeros+(1:opts.n_carriers);
    opts.type = mkms_type.full;
end

