% CIM_CC_MEX_SCRIPT   Generate MEX-function cim_cc_mex from cim_cc.
% 
% Script generated from project 'cim_cc.prj' on 10-May-2019.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.SaturateOnIntegerOverflow = false;
cfg.StackUsageMax = 900000;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'cim_cc'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg cim_cc -args ARGS{1} -nargout 2

