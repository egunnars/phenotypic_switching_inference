%% Background information

% This is a script which generates parameter regimes and artificial data as
% described in Appendix "Generation of artificial data" of the accompanying
% paper. It then implements the statistical model from Section "Estimation
% for cell number data" of the paper with E_(i,l)^(num) = 0, i.e. no
% measurement error.

% In the script, the parameters of the branching process model from Section
% "Multitype branching process" of the accompanying paper are arranged in a
% K x (K+1) matrix as follows:

% 1. The first column holds the birth rates of the phenotypes, i.e. entry
% (k,1) is the birth rate b_k of the k-th phenotype.

% 2. The second column holds the net birth rates of the phenotypes, i.e.
% entry (k,2) is the net birth rate lambda_k of the k-th phenotype.

% 3. The remaining columns hold the switching rates. The k-th row from the
% third column and onwards holds the switching rates from type-k to type-m
% in ascending order of m.

% As an example, consider the following parameter matrix for a three-type
% model:

% [0.9 0.5 0.02 0.04]
% [0.4 0.2 0.05 0.01]
% [0.6 0.3 0.03 0.08]

% Here, b_1 = 0.9, b_2 = 0.4 and b_3 = 0.6 are the birth rates of the
% phenotypes, lambda_1 = 0.5, lambda_2 = 0.2 and lambda_3 = 0.3 are the net
% birth rates of the phenotypes, nu_(12) = 0.02 and nu_(13) = 0.04 are the
% switching rates from type-1 to the other types, nu_(21) = 0.05 and
% nu_(23) = 0.01 are the switching rates from type-2, and nu_(31) = 0.03
% and nu_(32) = 0.08 are the switching rates from type-3.

%% Output

% The output of the script is the following:

% params_true_num is a K x (K+1) x nregimes tensor, which holds the true
% parameters for the artificial data. The true parameters are generated as
% described in Appendix "Generation of artifical data" of the accompanying
% paper. Arrangement of parameters is as described above.

% data_num is an I x K x L x R x nregimes x niter tensor which holds
% the artificial data.

% params_simple_num is a K x (K+1) x nregimes x niter tensor, which holds
% the estimates for the model parameters under the deterministic population
% model presented in Appendix "Implementation in MATLAB" of the
% accompanying paper.
 
% params_mle_num is a K x (K+1) x nregimes x niter tensor, which holds the
% estimates for the model parameters under the statistical model from
% Section "Estimation for cell number data", with no measurement error term.

% cvs_num is a K x (K+1) x nregimes x 2 tensor. cvs_num(i,j,n,1) holds the
% coefficient of variation (CV) for the simple deterministic estimates of
% parameter (i,j) under parameter regime n. cvs_num(i,j,n,2) holds the CV
% for the estimates of parameter (i,j) under the statistical model from
% Section "Estimation for cell number data".

%% Input parameters

% A random seed can be fixed for reproducibility.
rndseed = 1; rng(rndseed);

% K is the number of phenotypes in the model.
K = 2;
% I is the number of initial conditions used in the cell line experiments.
I = 2;
% R is the number of experimental replicates performed.
R = 3;
% T is the row vector of experimental timepoints.
T = [1:1:6];
% L is the number of experimental timepoints.
L = size(T,2);

% nregimes is the number of parameter regimes to generate
nregimes = 2;
% niter is the number of artfifical datasets to generate for each parameter
% regime.
niter = 10;

% If time_lapse = 0, the data is generated as endpoint data. If timelapse =
% 1, the data is generated as sequential data. See Section "Expmeriments
% and data" of the accompanying paper.
time_lapse = 0;
% If dead_option = 1, the number of dead cells at each time point is used
% in the parameter estimation. See Section "Improving identifiability of
% the rates of cell division and cell death" of the accompanying paper.
dead_option = 0;

% C_num is a cell array, where C_num{i} is a K x J_i matrix which
% allows the user to reduce the experimental data under the i-th initial
% condition. This option can be useful for reducible switching dynamics.
% See Appendices "Estimation for reducible switching dynamics" and
% "Implementation in MATLAB" of the accompanying paper.
C_num = cell(I);
for i=1:I
    C_num{i} = eye(K,K);
end

% lb_num is a K x (K+1) matrix which holds lower bounds for the branching
% process model parameters. The arrangement of parameters is as described
% at the top of the script. For example, the first column holds lower
% bounds for the birth rates of the phenotypes.
lb_num = [0,-20,0;0,-20,0];
% ub_num is a K x (K+1) matrix which holds upper bounds for the branching
% process model parameters.
ub_num = [10,20,10;10,20,10];
% Aineq_num is a K x (K+1) x n tensor, where n is the number of linear
% inequality constraints on the branching process model parameters imposed
% by the user. Each inequality constraint is of the form
% a_1x_1 + ... + a_Mx_M <= b,
% where x_1,...,x_M are the branching process model parameters.
% For each m = 1,...,n, Aineq_num(:,:,m) holds the coefficients of the
% left-hand side of the m-th inequality constraint. The arrangement of
% parameters is as described at the top of the script.
Aineq_num = [];
% bineq_num is an n x 1 vector, where n is the number of linear inequality
% constraints. For each m = 1,...,n, bineq_num(m) holds the coefficient of
% the right-hand side of the m-th inequality constraint.
bineq_num = [];
% Aeq_num is a K x (K+1) x n tensor, where n is the number of linear
% equality constraints on the branching process model parameters imposed
% by the user. Each equality constraint is of the form
% a_1x_1 + ... + a_Mx_M = b,
% where x_1,...,x_M are the branching process model parameters.
% For each m = 1,...,n, Aeq_num(:,:,m) holds the coefficients of the
% left-hand side of the m-th equality constraint.
Aeq_num = [];
% beq_num is an n x 1 vector, where n is the number of linear equality
% constraints. For each m = 1,...,n, beq_num(m) holds the coefficient of
% the right-hand side of the m-th equality constraint.
beq_num = [];
% In this script, neither an inequality constraint nor an equality
% constraint is implemented.

% By default, the MLE optimization problem is solved once, using an initial
% guess based on simple deterministic parameter estimates. Here, the user
% is given the option to request that the optimization problem is solved
% for multiple initial guesses, each based on the simple deterministic
% estimates. See Appendix "Implementation in MATLAB" of the accompanying
% paper.
nopt_simple_mle = 0;
% Here, the user is given the option to request that the MLE optimization
% problem for is solved for multiple initial guesses, where the initial
% guesses are chosen in a random fashion as described in Appendix
% "Generation of artificial data" of the accompanying paper.
nopt_random_mle = 0;

%% Setup of data structures

data_num_proj = zeros(I,K,L,R);
dead_proj = zeros(I,1,L,R);

params_simple_num = zeros(K,K+1,nregimes,niter);
params_mle_num = zeros(K,K+1,nregimes,niter);
scale_mle_num = zeros(K,K+1,nregimes,niter);
optlikelihood_num = zeros(nregimes,niter);
optlikelihood_simple_num_2 = zeros(nregimes,niter);
feasible_num = zeros(nregimes,niter);

means_num = zeros(K,K+1,nregimes,2);
cvs_num = zeros(K,K+1,nregimes,2);
means_num_rel = zeros(K,K+1,nregimes);
cvs_num_rel = zeros(K,K+1,nregimes);

runs_num = zeros(nregimes,niter);

%% Generation of model parameters and initial conditions for experiments
[params_true, params_true_num, params_true_frac, scale_true_num, scale_true_frac, N, Ntot, alpha] = params_generation(K,nregimes);

%% Generation of synthetic data via simulation
tic
[data_num, data_frac, dead] = data_generation(T,N,K,I,L,R,niter,params_true_num,time_lapse);
time_simul = toc;
time_simul = time_simul/(nregimes*niter);

if dead_option == 0
   dead = []; 
end

%% Estimation

tol = 1e-3;

tic
for reg=1:nregimes
    for k=1:niter

        data_num_proj = data_num(:,:,:,:,reg,k);
        if isempty(dead)
            dead_proj = [];
        else
            dead_proj = dead(:,:,:,:,reg,k);
        end
        
        end_cond = 0;       
        while end_cond == 0
            params_simple_num(:,:,reg,k) = estimation_simple_num(data_num_proj,N(:,:,reg),T);   
            [x0_num, scale_num] = initialguess_num(dead_proj,N(:,:,reg),T,params_simple_num(:,:,reg,k),lb_num,ub_num);
            [params_mle_num(:,:,reg,k), scale_mle_num(:,:,reg,k), optlikelihood_num(reg,k), feasible_num(reg,k)] ...
             = estimation_mle_num(data_num_proj,dead_proj,C_num,N(:,:,reg),T,x0_num,scale_num,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,tol);
            optlikelihood_simple_num_2(reg,k) = neg_log_likelihood_num(x0_num,ones(K,K+1),C_num,N(:,:,reg),data_num_proj,dead_proj,T);
            runs_num(reg,k) = runs_num(reg,k)+1;
            
            if sum(sum(abs(params_simple_num(:,2:end,reg,k)-params_mle_num(:,2:end,reg,k))./min(abs(params_simple_num(:,2:end,reg,k)),abs(params_mle_num(:,2:end,reg,k)))>1)) == 0 ...
                && optlikelihood_num(reg,k) <= optlikelihood_simple_num_2(reg,k)
                end_cond = 1;
            end
        end
            
        if nopt_simple_mle > 0 || nopt_random_mle > 0
            for n = 1:nopt_simple_mle
                [x0_num, scale_num] = initialguess_num(dead,N,T,params_simple_num,lb_num,ub_num);
                [params_mle_temp, scale_temp, optvalue_temp, feasible_temp] ...
                    = estimation_mle_num(data_num,dead,C_num,N,T,x0_num,scale_num,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,tol);
            
                if optvalue_temp < optlikelihood_num
                    params_mle_num = params_mle_temp;
                    scale_mle_num = scale_temp;
                    optlikelihood_num = optvalue_temp;
                    feasible_num = feasible_temp;
                end     
            end
            for n = 1:nopt_random_mle
                [x0_num, scale_num] = initialguess_num(dead,N,T,[],lb_num,ub_num);
                [params_mle_temp, scale_temp, optvalue_temp, feasible_temp] ...
                    = estimation_mle_num(data_num,dead,C_num,N,T,x0_num,scale_num,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,tol);
            
                if optvalue_temp < optlikelihood_num
                    params_mle_num = params_mle_temp;
                    scale_mle_num = scale_temp;
                    optlikelihood_num = optvalue_temp;
                    feasible_num = feasible_temp;
                end     
            end
        end
    end
end
time_mle_num = toc;
time_mle_num = time_mle_num/(nregimes*niter);

%% Summary statistics

%Rearrange parameters for comparison with cell fraction estimation
params_true_num_rel = params_true_num;
params_mle_num_rel = params_mle_num;
for reg=1:nregimes
   params_true_num_rel(2:end,2,reg) =  params_true_num(2:end,2,reg)-params_true_num(1,2,reg);
   for k=1:niter
       params_mle_num_rel(2:end,2,reg,k) = params_mle_num(2:end,2,reg,k)-params_mle_num(1,2,reg,k);
   end
end
    
%Computation of summary statistics
for reg=1:nregimes
    for i=1:K
        for j=1:K+1
            means_num(i,j,reg,1) = abs(mean(params_simple_num(i,j,reg,:))/params_true_num(i,j,reg)-1);
            means_num(i,j,reg,2) = abs(mean(params_mle_num(i,j,reg,:))/params_true_num(i,j,reg)-1);
            means_num_rel(i,j,reg) = abs(mean(params_mle_num_rel(i,j,reg,:))/params_true_num_rel(i,j,reg)-1);
            cvs_num(i,j,reg,1) = abs(std(params_simple_num(i,j,reg,:))/mean(params_simple_num(i,j,reg,:)));
            cvs_num(i,j,reg,2) = abs(std(params_mle_num(i,j,reg,:))/mean(params_mle_num(i,j,reg,:)));
            cvs_num_rel(i,j,reg) = abs(std(params_mle_num_rel(i,j,reg,:))/mean(params_mle_num_rel(i,j,reg,:)));
        end
    end
end
