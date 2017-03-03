%                               FASTA.M
%      This method solves the problem
%                        minimize f(Ax)+g(x)
%   Where A is a matrix, f is differentiable, and both f and g are convex.  
%   The algorithm is an adaptive/accelerated forward-backward splitting.  
%   The user supplies function handles that evaluate 'f' and 'g'.  The user 
%   also supplies a function that evaluates the gradient of 'f' and the
%   proximal operator of 'g', which is given by
%                proxg(z,t) = argmin t*g(x)+.5||x-z||^2.
%
%  Inputs:
%    A     : A matrix (or optionally a function handle to a method) that 
%             returns A*x
%    At    : The adjoint (transpose) of 'A.' Optionally, a function handle
%             may be passed.  
%    gradf : A function of z, computes the gradient of f at z
%    proxg : A function of z and t, the proximal operator of g with
%             stepsize t.  
%    x0    : The initial guess, usually a vector of zeros
%    f     : A function of x, computes the value of f
%    g     : A function of x, computes the value of g
%    opts  : An optional struct with options.  The commonly used fields
%             of 'opts' are:
%               maxIters : (integer, default=1e4) The maximum number of iterations
%                               allowed before termination.
%               tol      : (double, default=1e-3) The stopping tolerance.  
%                               A smaller value of 'tol' results in more
%                               iterations.
%               verbose  : (boolean, default=false)  If true, print out
%                               convergence information on each iteration.
%               recordObjective:  (boolean, default=false) Compute and
%                               record the objective of each iterate.
%               recordIterates :  (boolean, default=false) Record every
%                               iterate in a cell array.
%            To use these options, set the corresponding field in 'opts'. 
%            For example:  
%                      >> opts.tol=1e-8;
%                      >> opts.maxIters = 100;
%
%  Outputs:
%    sol  : The approximate solution
%    outs : A struct with convergence information
%    opts : A complete struct of options, containing all the values 
%           (including defaults) that were used by the solver.
%
%   For more details, see the FASTA user guide, or the paper "A field guide
%   to forward-backward splitting with a FASTA implementation."
%
%   Copyright: Tom Goldstein, 2014.


function [ sol, outs, opts ] = fasta(A, At, f, gradf, g, proxg, x0, opts )

 %%  Check whether we have function handles or matrices
 if ~isnumeric(A)
     assert(~isnumeric(At),'If A is a function handle, then At must be a handle as well.')
 end
 %  If we have matrices, create functions so we only have to treat one case
 if isnumeric(A)
     At = @(x)A'*x;
     A = @(x) A*x;
 end
 
 %% Check preconditions, fill missing optional entries in 'opts'
 if ~exist('opts','var') % if user didn't pass this arg, then create it
     opts = [];
 end
 opts = setDefaults(opts,A,At,x0,gradf); % fill default values for options
 % Verify that At=A'
 checkAdjoint(A,At,x0);
 
 %if opts.verbose
 %   fprintf('%sFASTA:\tmode = %s\n\tmaxIters = %i,\ttol = %1.2d\n',...
 %                      opts.stringHeader, opts.mode, opts.maxIters, opts.tol);
 %end
 
 %% Record some frequently used information from opts
 tau1      = opts.tau;                  % initial stepsize
 max_iters = opts.maxIters;             % maximum iterations before automatic termination
 W         = opts.window;               % lookback window for non-montone line search
 
 %% Allocate memory
 residual   = zeros(max_iters,1);       %  Residuals
 normalizedResid = zeros(max_iters,1);  %  Normalized residuals
 taus       = zeros(max_iters,1);       %  Stepsizes
 fVals      = zeros(max_iters,1);       %  The value of 'f', the smooth objective term
 objective  = zeros(max_iters+1,1);     %  The value of the objective function (f+g)
 funcValues = zeros(max_iters,1);       %  Values of the optional 'function' argument in 'opts'
 totalBacktracks = 0;                   %  How many times was backtracking activated?
 backtrackCount  = 0;                   %  Backtracks on this iterations

 %% Intialize array values
 x1       = x0;
 d1       = A(x1);
 f1       = f(d1);
 fVals(1) = f1;
 gradf1   = At(gradf(d1));
 
 %%  Initialize additional storage required for FISTA
 if opts.accelerate
     x_accel1 = x0;
     d_accel1 = d1;
     alpha1   = 1;
 end
 
 %  To handle non-monotonicity
 maxResidual       = -Inf; %  Stores the maximum value of the residual that has been seen. Used to evaluate stopping conditions.
 minObjectiveValue =  Inf;  %  Stores the best objective value that has been seen.  Used to return best iterate, rather than last iterate
 
 %  If user has chosen to record objective, then record initial value
 if opts.recordObjective %  record function values
     objective(1) = f1+g(x0); 
 end
 
 tic;   % Begin recording solve time
 %% Begin Loop
 for i = 1:max_iters
     %%  Rename iterates relative to loop index.  "0" denotes index i, and "1" denotes index i+1
     x0=x1;              % x_i <--- x_{i+1}
     gradf0 = gradf1;    % gradf0 is now $\nabla f (x_i)$ 
     tau0 = tau1;        % \tau_i <--- \tau_{i+1}

     %%  FBS step: obtain x_{i+1} from x_i
     x1hat = x0 - tau0*gradf0; % Define \hat x_{i+1}
     x1 = proxg(x1hat,tau0);   % Define x_{i+1}
    
     %%  Non-monotone backtracking line search
     Dx = x1-x0;
     d1 = A(x1); 
     f1 = f(d1);
     if opts.backtrack
         M = max( fVals(max(i-W,1):max(i-1,1)) );  % Get largest of last 10 values of 'f'
         backtrackCount=0;
         %  Note: 1e-12 is to quench rounding errors
         while f1-1e-12> M+real(dot(Dx(:),gradf0(:)))+norm(Dx(:))^2/(2*tau0) && backtrackCount<20  % the backtracking loop
          tau0 = tau0*opts.stepsizeShrink;    % shrink stepsize
          x1hat = x0 - tau0*gradf0; % redo the FBS
          x1 = proxg(x1hat,tau0);
          d1 = A(x1);
          f1 = f(d1);
          Dx = x1-x0; 
          backtrackCount = backtrackCount+1;
         end
         totalBacktracks = totalBacktracks+backtrackCount;
     end
     
    % if opts.verbose && backtrackCount>10
    %       fprintf('%s\tWARNING: excessive backtracking (%d steps), current stepsize is %f\n',...
     %                  opts.stringHeader, backtrackCount, tau0);
    % end
     
     %% Record convergence information
     taus(i) = tau0; % stepsize
     residual(i)  =  norm(Dx(:))/tau0; % Estimate of the gradient, should be zero at solution
     maxResidual = max(maxResidual, residual(i));
     normalizer = max(norm(gradf0(:)),norm(x1(:)-x1hat(:))/tau0)+opts.eps_n;
     normalizedResid(i) = residual(i)/normalizer;  % Normalized residual:  size of discrepancy between the two derivative terms, divided by the size of the terms
     fVals(i) = f1;
     funcValues(i) = opts.function(x0);
     if opts.recordObjective   %  Record function values
        objective(i+1) = f1+g(x1); 
        newObjectiveValue = objective(i+1);
     else
        newObjectiveValue = residual(i); %  Use the residual to evalue quality of iterate if we don't have objective
     end
     
     if opts.recordIterates   %  Record function values
        iterates{i} = x1;
     end
     
     if newObjectiveValue<minObjectiveValue  % Methods is non-monotone:  Make sure to record best solution
        bestObjectiveIterate = x1;
        minObjectiveValue = newObjectiveValue; 
     end
     
   %  if opts.verbose>1
    %  fprintf('%s%d: resid = %0.2d, backtrack = %d, tau = %d', opts.stringHeader, i, residual(i), backtrackCount,tau0);
    %      if opts.recordObjective
    %           fprintf(', objective = %d\n', objective(i+1));
    %      else
     %          fprintf('\n');
     %     end
    %  end


    %% Test stopping criteria
    %  If we stop, then record information in the output struct
     if opts.stopNow(x1,i,residual(i),normalizedResid(i),maxResidual,opts) || i>=max_iters
         outs = [];
         outs.solveTime = toc;
         outs.residuals = residual(1:i);
         outs.stepsizes = taus(1:i);
         outs.normalizedResiduals = normalizedResid(1:i);
         outs.objective = objective(1:i);
         outs.funcValues = funcValues(1:i);
         outs.backtracks = totalBacktracks;
         outs.L = opts.L;
         outs.initialStepsize = opts.tau;
         outs.iterationCount = i;
         if ~opts.recordObjective
            outs.objective = 'Not Recorded';
         end
         if opts.recordIterates
            outs.iterates = iterates;
         end
         sol = bestObjectiveIterate;
   %      if opts.verbose
   %          fprintf('%s\tDone:  time = %0.3f secs, iterations = %i\n',opts.stringHeader, toc, outs.iterationCount);
   %      end
         return;
     end
     
     if opts.adaptive && ~opts.accelerate
         %% Compute stepsize needed for next iteration using BB/spectral method
         gradf1 = At(gradf(d1));
         Dg = gradf1+(x1hat-x0)/tau0;% Delta_g, note that Delta_x was recorded above during backtracking
         dotprod = real(dot(Dx(:),Dg(:)));
         tau_s = norm(Dx(:))^2/ dotprod;  %  First BB stepsize rule
         tau_m = dotprod / norm(Dg(:))^2; %  Alternate BB stepsize rule      
         tau_m = max(tau_m,0);
         if 2*tau_m > tau_s   %  Use "Adaptive"  combination of tau_s and tau_m
             tau1 = tau_m;
         else
             tau1 = tau_s - .5*tau_m;  %  Experiment with this param
         end 
         if tau1 <=0 || isinf(tau1) || isnan(tau1)      %  Make sure step is non-negative
             tau1 = tau0*1.5;  % let tau grow, backtracking will kick in if stepsize is too big
         end
     end

     if opts.accelerate
         %% Use FISTA-type acceleration
         x_accel0 = x_accel1;  %  Store the old iterates
         d_accel0 = d_accel1;
         alpha0=alpha1; 
         x_accel1 = x1;
         d_accel1 = d1;
         %  Check to see if the acceleration needs to be restarted
         if opts.restart && (x0(:)-x1(:))'*(x1(:)-x_accel0(:))>0
             alpha0=1;
         end
         %  Calculate acceleration parameter
         alpha1 = (1+sqrt(1+4*alpha0^2))/2;
         %  Over-relax/predict
         x1 = x_accel1+(alpha0-1)/alpha1*(x_accel1-x_accel0);
         d1 = d_accel1+(alpha0-1)/alpha1*(d_accel1-d_accel0);
         %  Compute the gradient needed on the next iteration
         gradf1 = At(gradf(d1));
         fVals(i) = f(d1);
         tau1 = tau0;
     end
     
     if ~opts.adaptive && ~opts.accelerate
          gradf1 = At(gradf(d1));
          tau1 = tau0;
     end
   
 end


return
function checkAdjoint(A,At,x)
 x = randn(size(x));
 Ax = A(x);
 y = randn(size(Ax));
 Aty = At(y);
 
 innerProduct1 = Ax(:)'*y(:);
 innerProduct2 = x(:)'*Aty(:);
 error = abs(innerProduct1-innerProduct2)...
     /max(abs(innerProduct1),abs(innerProduct2));
 assert(error<1e-9,'"At" is not the adjoint of "A".  Check the definitions of these operators.'); 
return

