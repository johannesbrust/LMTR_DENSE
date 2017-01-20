% This is a script to minimize functions from the CUTEst test set using 
% the Eig(inf,2) with a multiple of identity or denso initial matrix.
% This script is for one function only.
%
% Last modified: 11/17/16, J.B
%
% For details, see:
% O.Burdakov, L. Gong, Y. Yuan, S. Zikrin, On Efficiently Combining Limited
% Memory and Trust-Region Techniques, technical report LiTH-MAT-R-2013/13-SE,
% Department of Mathematics, Link�ping University, 2013.
% http://liu.diva-portal.org/smash/record.jsf?pid=diva2%3A667359

clc;
clear;

parentpath = [cd(cd('..')) '/main'];
addpath(parentpath)

% Set input parameters:
params=struct;
params.m = 5;  % number of L-BFGS updates
params.gtol = 1e-5;  % exit if ||g||_2<gtol*max(1,||x||_2)
params.ranktol = 1e-7;  % tolerance for establishing rank of V
params.dflag = 0;  % display parameter, 1 if to display information per iteration
params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
params.ftol=1e-11;  % tolerance on relative function reduction

params.gamma_opt=1;  % how to compute gamma_perp for the dense initial matrix.
params.gamma_store=1;  % how to compute gamma_perp for the dense initial matrix.

% Set maximum iterations, J.B 20/10
params.maxit = 10000;

CUTEst_init  % initialize CUTEr, see appropriate documentation 
%fid =fopen('cutest_list.txt','r');  % read file that contains CUTEr test problems

problem = 'BDQRTIC';


eval(['!runcutest -p matlab -D ' problem]);
prob    = cutest_setup();
x0      = prob.x; 

% Dense initial matrix
tic;
[x1,f1,out1,gamms1] = LMTR_EIG_inf_2_DENSE(@cutest_fun,x0,params);
t1 = toc;

% Multiple of identity
tic;
[x2,f2,out2] = LMTR_EIG_inf_2(@cutest_fun,x0,params);
t2 = toc;

norm(x1-x2)
norm(f1-f2)

