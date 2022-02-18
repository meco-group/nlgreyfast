function phases = mkfullms_random(opts)
    % Generate random phases for multisine (generate the signal itself using mkfullms).
    phases = 2*pi*(rand(opts.n_carriers,1)-0.5);
