% KRASKOVMI_CC_MEX_SCRIPT   Generate MEX-function KraskovMI_cc_mex from
%  KraskovMI_cc.
% 
% Script generated from project 'KraskovMI_cc.prj' on 10-May-2019.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'KraskovMI_cc'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{3} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg KraskovMI_cc -args ARGS{1} -nargout 1

