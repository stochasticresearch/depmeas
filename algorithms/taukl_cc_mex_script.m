% TAUKL_CC_MEX_SCRIPT   Generate MEX-function taukl_cc_mex from taukl_cc.
% 
% Script generated from project 'taukl_cc.prj' on 10-May-2019.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'taukl_cc'.
ARGS = cell(1,1);
ARGS{1} = cell(5,1);
ARGS{1}{1} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg taukl_cc -args ARGS{1} -nargout 1

