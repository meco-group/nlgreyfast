function solver = sysid_gauss_newton(e,nlp,V,fopts)
  % Helper function to create nlpsol object
  if nargin < 4, fopts = struct; end
  if(~isfield(fopts,"generate_on")) fopts.generate_on = false;  end
  if(~isfield(fopts,"compiler_flags")) fopts.compiler_flags = ''; end
  if(~isfield(fopts,"ipopt")) fopts.ipopt = struct; end
  if(~isfield(fopts,"hessian")) fopts.hessian = 'gaussnewton'; end
  
  J = jacobian(e,V);
  H = triu(J'*J);
  sigma = casadi.MX.sym('sigma');
  
  io = struct;
  io.x = V;
  io.lam_f = sigma;
  io.hess_gamma_x_x = sigma*H;
  
  opts = struct;
  opts.jit = false;
  %opts.dump_in = true;
  %opts.dump = true;
  
  if fopts.compiler_flags
    opts.jit = true;
    opts.compiler='shell';
    opts.jit_options.verbose = true;
    opts.jit_options.compiler_flags = compiler_flags;
    opts.jit_options.compiler = 'ccache gcc';
    opts.jit_cleanup = false;
    opts.jit_temp_suffix = false;
  else
    opts.jit = false;
  end
  
  if strcmp(fopts.hessian,'gaussnewton')
    hessLag = casadi.Function('nlp_hess_l',io,{'x','p','lam_f','lam_g'}, {'hess_gamma_x_x'},opts);
    opts.hess_lag = hessLag;
  elseif strcmp(fopts.hessian,'exact')
    %if it's exact then there's nothing to do -- if opts.hess_lag is not defined, it'll be like that.
  else
    error('unknown value for fopts.hessian');
  end
  opts.ipopt = fopts.ipopt;
  %opts.ipopt.bound_push=0.1;
  %opts.ipopt.max_cpu_time = 5*60;
  %opts.ipopt.tol = 1e-16;
  %opts.ipopt.max_iter = 500;
  %cbarg = casadi.MX.sym('cbarg');
  %opts.iteration_callback = casadi.Function('devnull',{cbarg},{0});
  solver = casadi.nlpsol('solver','ipopt', nlp, opts);
  if fopts.generate_on
      casadi.Function('generated_f',{nlp.x},{nlp.f},{'x'},{'obj'}).generate('nlidss_f.c',struct('with_header',true))
      casadi.Function('generated_j',{nlp.x},{J},{'x'},{'j'}).generate('nlidss_j.c',struct('with_header',true))
      casadi.Function('generated_jt',{nlp.x},{J.'},{'x'},{'jt'}).generate('nlidss_jt.c',struct('with_header',true))
      casadi.Function('generated_e',{nlp.x},{e},{'x'},{'e'}).generate('nlidss_e.c',struct('with_header',true))
  end
  
end
