%  This script shows how to use FASTA to solve the Lasso problem:
%                min  .5||Ax-b||^2 s.t. |x|<mu
%  Where A is an MxN matrix, b is an Mx1 vector of measurements, and x is
%  the Nx1 vector of unknowns.  The parameter 'mu' controls the strength of
%  the regularizer.


function [x] = fastiter(A,b,mu,x0)
%%  OPTIONAL:  give some extra instructions to FASTA using the 'opts' struct
opts = [];
%opts.tol = 1e-8;  % Use super strict tolerance
opts.recordObjective = true; %  Record the objective function so we can plot it
opts.verbose = true;
opts.stringHeader='    ';      % Append a tab to all text output from FISTA.  This option makes formatting look a bit nicer. 


%%  Call the solver 3 times
% Default behavior: adaptive stepsizes
[sol, outs_adapt] = fasta_sparseLeastSquares(A,A',b,mu,x0, opts);
x = sol;
opts.accelerate = false;
opts.adaptive = false;
[sol1, outs_fbs] = fasta_sparseLeastSquares(A,A',b,mu,x0, opts);
plotresidual;