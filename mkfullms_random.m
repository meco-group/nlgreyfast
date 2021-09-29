function phases = mkfullms_random(opts)
    phases = 2*pi*(rand(opts.n_carriers,1)-0.5);
