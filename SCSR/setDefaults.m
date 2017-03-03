%% Fill in the struct of options with the default values
function opts = setDefaults(opts,A,At,x0,gradf)


%  maxIters: The maximum number of iterations
if ~isfield(opts,'maxIters')
    opts.maxIters = 30;
end

% tol:  The relative decrease in the residuals before the method stops
if ~isfield(opts,'tol') % Stopping tolerance
    opts.tol = 1e-3;
end

% verbose:  If 'true' then print status information on every iteration
if ~isfield(opts,'verbose')   
    opts.verbose = false;
end

% recordObjective:  If 'true' then evaluate objective at every iteration
if ~isfield(opts,'recordObjective')   
    opts.recordObjective = false;
end

% recordIterates:  If 'true' then record iterates in cell array
if ~isfield(opts,'recordIterates')   
    opts.recordIterates = false;
end

% adaptive:  If 'true' then use adaptive method.
if ~isfield(opts,'adaptive')    %  is Adaptive?
    opts.adaptive = false;
end

% accelerate:  If 'true' then use FISTA-type adaptive method.
if ~isfield(opts,'accelerate')    %  is Accelerated?
    opts.accelerate = true;
end

% restart:  If 'true' then restart the acceleration of FISTA.
%   This only has an effect when opts.accelerate=true
if ~isfield(opts,'restart')    %  use restart?
    opts.restart = true;
end

% backtrack:  If 'true' then use backtracking line search
if ~isfield(opts,'backtrack')    
    opts.backtrack = false;
end

% stepsizeShrink:  Coefficient used to shrink stepsize when backtracking
% kicks in
if ~isfield(opts,'stepsizeShrink')    
    opts.stepsizeShrink = 0.2;      % The adaptive method can expand the stepsize, so we choose an aggressive value here
    if ~opts.adaptive || opts.accelerate
         opts.stepsizeShrink = 0.5; % If the stepsize is monotonically decreasing, we don't want to make it smaller than we need
    end
end

%  Create a mode string that describes which variant of the method is used
opts.mode = 'plain';
if opts.adaptive
     opts.mode = 'adaptive';
end
if opts.accelerate
    if opts.restart
        opts.mode = 'accelerated(FISTA)+restart';
    else
        opts.mode = 'accelerated(FISTA)';
    end
end


% W:  The window to look back when evaluating the max for the line search
if ~isfield(opts,'window') % Stopping tolerance
    opts.window = 10;
end

% eps_r:  Epsilon to prevent ratio residual from dividing by zero
if ~isfield(opts,'eps_r') % Stopping tolerance
    opts.eps_r = 1e-8;
end

% eps_n:  Epsilon to prevent normalized residual from dividing by zero
if ~isfield(opts,'eps_n') % Stopping tolerance
    opts.eps_n = 1e-8;
end

%  L:  Lipschitz constant for smooth term.  Only needed if tau has not been
%   set, in which case we need to approximate L so that tau can be
%   computed.
if (~isfield(opts,'L') || opts.L<=0) && (~isfield(opts,'tau') || opts.tau<=0)
    x1 = randn(size(x0));
    x2 = randn(size(x0));
    gradf1 = At(gradf(A(x1)));
    gradf2 = At(gradf(A(x2)));
    opts.L = norm(gradf1(:)-gradf2(:))/norm(x2(:)-x1(:));
    opts.L = max(opts.L,1e-6);
    opts.tau = 2/opts.L/10;
end
assert(opts.tau>0,['Invalid step size: ' num2str(opts.tau)]);

%  Set tau if L was set by user
if(~isfield(opts,'tau') || opts.tau<=0)
    opts.tau = 1.0/opts.L;
else
    opts.L = 1/opts.tau;
end

% function:  An optional function that is computed and stored after every
% iteration
if ~isfield(opts,'function')          % This functions gets evaluated on each iterations, and results are stored
    opts.function = @(x) 0;
end

% stringHeader:  Append this string to beginning of all output
if ~isfield(opts,'stringHeader')          % This functions gets evaluated on each iterations, and results are stored
    opts.stringHeader = '';
end

%  The code below is for stopping rules
%  The field 'stopNow' is a function that returns 'true' if the iteration
%  should be terminated.  The field 'stopRule' is a string that allows the
%  user to easily choose default values for 'stopNow'.  The default
%  stopping rule terminates when the relative residual gets small.
if isfield(opts,'stopNow') 
    opts.stopRule = 'custom';
end

if ~isfield(opts,'stopRule') 
    opts.stopRule = 'hybridResidual';
end

if strcmp(opts.stopRule,'residual')
    opts.stopNow = @(x1,iter,resid,normResid,maxResidual,opts) resid < opts.tol; 
end

if strcmp(opts.stopRule,'iterations')
    opts.stopNow = @(x1,iter,resid,normResid,maxResidual,opts) iter > opts.maxIters; 
end

% Stop when normalized residual is small
if strcmp(opts.stopRule,'normalizedResidual')
    opts.stopNow = @(x1,iter,resid,normResid,maxResidual,opts) normResid < opts.tol; 
end

% Divide by residual at iteration k by maximum residual over all iterations.
% Terminate when this ratio gets small.
if strcmp(opts.stopRule,'ratioResidual')
    opts.stopNow = @(x1,iter,resid,normResid,maxResidual,opts)   resid/(maxResidual+opts.eps_r) < opts.tol; 
end

% Default behavior:  Stop if EITHER normalized or ration residual is small
if strcmp(opts.stopRule,'hybridResidual')
    opts.stopNow = @(x1,iter,resid,normResid,maxResidual,opts) ...
        resid/(maxResidual+opts.eps_r) < opts.tol ...
        || normResid < opts.tol; 
end

assert(isfield(opts,'stopNow'),['Invalid choice for stopping rule: ' opts.stopRule ]);

return
