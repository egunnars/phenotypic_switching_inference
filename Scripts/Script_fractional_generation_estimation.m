%% Background information

% This is a script which generates parameter regimes and artificial data as
% described in Appendix "Generation of artificial data" of the accompanying
% paper. It then implements estimation under the statistical model from
% Section "Estimation for cell fraction data" of the accompanying paper
% with E_(i,l)^(frac) = 0, i.e. no measurement error.

% In the script, the parameters of the branching process model from Section
% "Multitype branching process" of the accompanying paper are arranged in a
% K x (K+1) matrix as follows:

% 1. The first column holds the death rates of the phenotypes, i.e. entry
% (k,1) is the death rate d_k of the k-th phenotype.

% 2. The first entry of the second column holds the net birth rate lambda_1
% for type-1 cells, while entry k for k >= 2 holds the relative net birth
% rate lambda_k-lambda_1.

% 3. The remaining columns hold the switching rates. The k-th row from the
% third column and onwards holds the switching rates from type-k to type-m
% in ascending order of m.

% As an example, consider the following parameter matrix for a three-type
% model:

% [1.0 0.6 0.02 0.04]
% [0.8 0.4 0.05 0.01]
% [0.7 0.3 0.03 0.08]

% Here, d_1 = 1.0, d_2 = 0.8 and d_3 = 0.7 are the death rates of the
% phenotypes, lambda_1 = 0.6 is the net birth rate of type-1 cells,
% lambda_2-lambda_1 = 0.4 and lambda_3-lambda_1 = 0.3 are the relative net
% birth rates rates of type-2 and type-3 cells, nu_(12) = 0.02 and nu_(13)
% = 0.04 are the switching rates from type-1 to the other types, nu_(21) =
% 0.05 and nu_(23) = 0.01 are the switching rates from type-2, and nu_(31)
% = 0.03 and nu_(32) = 0.08 are the switching rates from type-3.

%% Output

% The output of the script is the following:

% params_true_frac is a K x (K+1) x nregimes tensor, which holds the true
% parameters for the artificial data. The true parameters are generated as
% described in Appendix "Generation of artificial data" of the accompanying
% paper. Arrangement of parameters is as described above.

% data_frac is an I x K x L x R x nregimes x niter tensor which holds
% the artificial data.

% params_simple_frac is a K x (K+1) x nregimes x niter tensor, which holds
% simple deterministic estimates for the model parameters, obtained by
% solving the least squares problem presented in Appendix "Implementation
% in MATLAB" of the accompanying paper.
 
% params_mle_frac is a K x (K+1) x nregimes x niter tensor, which holds the
% estimates for the model parameters under the statistical model from
% Section "Estimation for cell fraction data" of the accompanying paper,
% with no measurement error term.

% cvs_frac is a K x (K+1) x nregimes x 2 tensor. cvs_frac(i,j,n,1) holds the
% coefficient of variation (CV) for the simple deterministic estimates of
% parameter (i,j) under parameter regime n. cvs_frac(i,j,n,2)
% holds the CV for the estimates of parameter (i,j) under the statistical
% model from Section "Estimation for cell fraction data".

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

% C_frac is a cell array, where C_frac{i} is a (K-1) x J_i matrix which
% allows the user to reduce the experimental data under the i-th initial
% condition. This option can be useful for reducible switching dynamics.
% See Appendices "Estimation for reducible switching dynamics" and
% "Implementation in MATLAB" of the accompanying paper.
C_frac = cell(I);
for i=1:I
    C_frac{i} = eye(K-1,K-1);
end

% lb_frac is a K x (K+1) matrix which holds lower bounds for the branching
% process model parameters. The arrangement of parameters is as described
% at the top of the script. For example, the first column holds lower
% bounds for the death rates of the phenotypes.
lb_frac = [0,-20,0;0,-20,0];
% ub_frac is a K x (K+1) matrix which holds upper bounds for the branching
% process model parameters.
ub_frac = [10,20,10;10,20,10];
% Aineq_frac is a K x (K+1) x n tensor, where n is the number of linear
% inequality constraints on the branching process model parameters imposed
% by the user. Each inequality constraints is of the form
% a_1x_1 + ... + a_Mx_M <= b,
% where x_1,...,x_M are the branching process model parameters.
% For each m = 1,...,n, Aineq_frac(:,:,m) holds the coefficients of the
% left-hand side of the m-th inequality constraint. The arrangement of
% parameters is as described at the top of the script.
Aineq_frac = [];
% bineq_frac is an n x 1 vector, where n is the number of linear inequality
% constraints. For each m = 1,...,n, bineq_frac(m) holds the coefficient of
% the right-hand side of the m-th inequality constraint.
bineq_frac = [];
% Aeq_frac is a K x (K+1) x n tensor, where n is the number of linear
% equality constraints on the branching process model parameters imposed
% by the user. Each equality constraints is of the form
% a_1x_1 + ... + a_Mx_M = b,
% where x_1,...,x_M are the branching process model parameters.
% For each m = 1,...,n, Aeq_frac(:,:,m) holds the coefficients of the
% left-hand side of the m-th equality constraint.
Aeq_frac = [];
% beq_frac is an n x 1 vector, where n is the number of linear equality
% constraints. For each m = 1,...,n, beq_frac(m) holds the coefficient of
% the right-hand side of the m-th equality constraint.
beq_frac = [];
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

data_frac_proj = zeros(I,K,L,R);

params_simple_frac = zeros(K,K+1,nregimes,niter);
params_mle_frac = zeros(K,K+1,nregimes,niter);
scale_mle_frac = zeros(K,K+1,nregimes,niter);
optlikelihood_simple_frac = zeros(nregimes,niter);
optlikelihood_frac = zeros(nregimes,niter);
optlikelihood_simple_frac_2 = zeros(nregimes,niter);
feasible_frac = zeros(nregimes,niter);

means_frac = zeros(K,K+1,nregimes,2);
cvs_frac = zeros(K,K+1,nregimes,2);
means_frac_rel = zeros(K,K+1,nregimes,1);
cvs_frac_rel = zeros(K,K+1,nregimes,1);

runs_frac = zeros(nregimes,niter);

%% Generation of model parameters and initial conditions for experiments
[params_true, params_true_num, params_true_frac, scale_true_num, scale_true_frac, N, Ntot, alpha] = params_generation(K,nregimes);

%% Generation of synthetic data via simulation
tic
[data, data_frac, dead] = data_generation(T,N,K,I,L,R,niter,params_true_num,time_lapse);
time_simul = toc;
time_simul = time_simul/(nregimes*niter);

%% Estimation

tol = 1e-3;

tic
for reg=1:nregimes
    for k=1:niter

        data_frac_proj = data_frac(:,:,:,:,reg,k);
        
        end_cond = 0;       
        while end_cond == 0
            [params_simple_frac(:,:,reg,k), optlikelihood_simple_frac(reg,k)] = estimation_simple_frac(data_frac_proj,alpha(:,:,reg),Ntot(1,:,reg),T,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac);
            [x0_frac, scale_frac] = initialguess_frac(params_simple_frac(:,:,reg,k),lb_frac,ub_frac);
            [params_mle_frac(:,:,reg,k), scale_mle_frac(:,:,reg,k), optlikelihood_frac(reg,k), feasible_frac(reg,k)] ...
             = estimation_mle_frac(data_frac_proj,C_frac,Ntot(1,:,reg),alpha(:,:,reg),T,x0_frac,scale_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);
            optlikelihood_simple_frac_2(reg,k) = neg_log_likelihood_frac(x0_frac,ones(K,K+1),C_frac,Ntot(1,:,reg),alpha(:,:,reg),data_frac_proj,T);
            runs_frac(reg,k) = runs_frac(reg,k)+1;
            
            if sum(sum(abs((params_simple_frac(:,3:end,reg,k)-params_mle_frac(:,3:end,reg,k)))./min(abs(params_simple_frac(:,3:end,reg,k)),abs(params_mle_frac(:,3:end,reg,k)))>5)) ...
                 + sum((abs((params_simple_frac(2:end,2,reg,k)-params_mle_frac(2:end,2,reg,k)))./min(abs(params_simple_frac(2:end,2,reg,k)),abs(params_mle_frac(2:end,2,reg,k)))>1)) <= 1 ...
                 && optlikelihood_frac(reg,k) <= optlikelihood_simple_frac_2(reg,k)
                    end_cond = 1;
            end
        end
        
        if nopt_simple_mle > 0 || nopt_random_mle > 0
            for n = 1:nopt_simple_mle
                [x0_frac, scale_frac] = initialguess_frac(params_simple_frac(:,:,reg,k),lb_frac,ub_frac);
                [params_mle_temp, scale_temp, optvalue_temp, feasible_temp] ...
                    = estimation_mle_frac(data_frac_proj,C_frac,Ntot(1,:,reg),alpha(:,:,reg),T,x0_frac,scale_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);
                runs_frac(reg,k) = runs_frac(reg,k)+1;

                if optvalue_temp < optlikelihood_frac(reg,k)
                    params_mle_frac(:,:,reg,k) = params_mle_temp;
                    scale_mle_frac(:,:,reg,k) = scale_temp;
                    optlikelihood_frac(reg,k) = optvalue_temp;
                    feasible_frac(reg,k) = feasible_temp;
                end     
            end
            for n = 1:nopt_random_mle
                [x0_frac, scale_frac] = initialguess_frac([],lb_frac,ub_frac);
                [params_mle_temp, scale_temp, optvalue_temp, feasible_temp] ...
                    = estimation_mle_frac(data_frac_proj,C_frac,Ntot(1,:,reg),alpha(:,:,reg),T,x0_frac,scale_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);
                runs_frac(reg,k) = runs_frac(reg,k)+1;

                if optvalue_temp < optlikelihood_frac(reg,k)
                    params_mle_frac(:,:,reg,k) = params_mle_temp;
                    scale_mle_frac(:,:,reg,k) = scale_temp;
                    optlikelihood_frac(reg,k) = optvalue_temp;
                    feasible_frac(reg,k) = feasible_temp;
                end     
            end
        end
    end
end
time_mle_frac = toc;
time_mle_frac = time_mle_frac/(nregimes*niter);
    
%% Summary statistics

%Rearrange parameters for comparison with cell number estimation
params_true_frac_rel = params_true_frac;
params_mle_frac_rel = params_mle_frac;
for reg=1:nregimes
   params_true_frac_rel(:,1,reg) = params_true(:,1,reg);
   params_true_frac_rel(2:end,2,reg) = params_true_frac(2:end,2,reg);
   for k=1:niter
       params_mle_frac_rel(1,1,reg,k) = params_mle_frac(1,2,reg,k)+params_mle_frac(1,1,reg,k);
       params_mle_frac_rel(2:end,1,reg,k) = params_mle_frac(2:end,1,reg,k)+params_mle_frac(2:end,2,reg,k)+params_mle_frac(1,2,reg,k);
   end
end

%Compute summary statistics
for reg=1:nregimes
    for i=1:K
        for j=1:K+1
            means_frac(i,j,reg,1) = abs(mean(params_simple_frac(i,j,reg,:))/params_true_frac(i,j,reg)-1);
            means_frac(i,j,reg,2) = abs(mean(params_mle_frac(i,j,reg,:))/params_true_frac(i,j,reg)-1);
            means_frac_rel(i,j,reg) = abs(mean(params_mle_frac_rel(i,j,reg,:))/params_true_frac_rel(i,j,reg)-1);
            cvs_frac(i,j,reg,1) = abs(std(params_simple_frac(i,j,reg,:))/mean(params_simple_frac(i,j,reg,:)));
            cvs_frac(i,j,reg,2) = abs(std(params_mle_frac(i,j,reg,:))/mean(params_mle_frac(i,j,reg,:)));
            cvs_frac_rel(i,j,reg) = abs(std(params_mle_frac_rel(i,j,reg,:))/mean(params_mle_frac_rel(i,j,reg,:)));
        end
    end
end
