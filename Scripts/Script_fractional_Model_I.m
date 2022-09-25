%% Background information

% This is a script for cell fraction data with measurement error.

% The script implements estimation under Model I of Section "Application:
% Transition between stem and non-stem cell states in SW620 colon cancer"
% of the accompanying paper on a publicly available cell fraction dataset
% from Yang et al. (2012). It produces numerical results displayed in
% Section "Application: Transition between stem and non-stem cell states in
% SW620 colon cancer."

% In the script, the parameters of the branching process model from Section
% "Multitype branching process model" of the accompanying paper are
% arranged in a K x (K+1) matrix as follows:

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

% params_simple_frac is a K x (K+1) matrix which holds simple deterministic
% estimates for the model parameters, obtained by solving the least squares
% problem presented in Appendix "Implementation in MATLAB" of the
% accompanying paper. Arrangement of parameters is as described above.

% params_mle_frac is a K x (K+1) matrix which holds the estimates for
% the model parameters under Model I.

% If confidence intervals are requested, ci_frac is a 2 x K x (K+1)
% tensor, where for each i = 1,...,K and j = 1,...,K+1, ci_frac(:,i,j)
% holds the lower and upper interval bounds for parameter (i,j), with the
% arrangement of parameters described at the top of the script.

% For example, ci_frac(:,1,3) will return the confidence interval for the
% switching rate nu_(1,2).

% If requested, the script also produces a visualization of how well the
% model fits the data. See further explanation below.

%% Input data

% A random seed can be fixed for reproducibility.
rndseed = 1; rng(rndseed);

% K is the number of phenotypes in the model.
K = 2;
% I is the number of initial conditions used in the cell line experiments.
I = 2;
% R is the number of experimental replicates performed.
R = 1;
% T is the row vector of experimental timepoints.
T = [2:2:24];
% L is the number of experimental timepoints.
L = size(T,2);

% data_frac is an I x K x L x R tensor that holds the experimental data.
% Entry number (i,k,l,r) is the fraction of cells of type-k in the r-th
% replicate of the experiment started by the i-th initial condition and
% ended at the l-th timepoint.
% Here, we have manually entered the data from Yang et al. (2012).
data_frac = zeros(I,K,L,R);
data_frac(1,1,:,1) = [0.8350,0.7868,0.7323,0.7150,0.6940,0.7055,0.7140,0.7135,0.7345,0.7123,0.6788,0.6567];
data_frac(1,2,:,1) = 1-data_frac(1,1,:,1);
data_frac(2,1,:,1) = [0.3090,0.4000,0.4740,0.5523,0.5365,0.5043,0.5483,0.5960,0.5978,0.5920,0.6273,0.6190];
data_frac(2,2,:,1) = 1-data_frac(2,1,:,1);

% Ntot is a 1 x I row vector contaning the total starting number of cells
% under each initial condition.
Ntot = [10^3,10^3];
% alpha is an I x K matrix where the i-th row contains the starting fraction
% of cells of each type under the i-th initial condition.
alpha = [1,0;0,1];
% N is an I x K matrix where the i-th row contains the starting number of
% cells of each type under the i-th initial condition.
N = Ntot'.*alpha;

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
lb_frac = [0,-0.5,0;0,-0.5,0];
% ub_frac is a K x (K+1) matrix which holds upper bounds for the branching
% process model parameters.
ub_frac = [1,0.5,0.5;1,0.5,0.5];
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

% x0_frac_def is a K x (K+1) x n tensor, which can be used to supply
% initial guesses for the MLE optimization problem. Here, n is the number
% of distinct initial guesses the user wishes to supply. The arrangement of
% parameters is as described at the top of the script.
x0_frac_def = [];

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

% ci_option is a K x (K+1) matrix holding zeros and ones. Setting
% ci_option(i,j) = 1 will compute a confidence interval for parameter (i,j)
% of the branching process model, where the arrangement of parameters is
% as described at the top of the script.
% In this case, a confidence interval is requested for the switching rates
% nu_(12) and nu_(13) and the relative net birth rate lambda_2-lambda_1.
ci_option = zeros(K,K+1);
ci_option(1,3) = 1;
ci_option(2,3) = 1;
ci_option(2,2) = 1;
% alpha_q sets the level for the confidence intervals.
alpha_q = 0.05;
% nopt_ci allows the user to request that confidence intervals are computed
% multiple times, starting from different initial guesses.
nopt_ci = 0;

% If data_vis_option is set to 1, plots are produced that show for each
% initial condition and each type how well the statistical model fits the
% data. More precisely, the mean prediction of the statistical model and
% 1-alpha_q confidence intervals, under the assumption that the MLE
% estimates are the true parameters, are compared with the data.
data_vis_option = 1;

%% Setup of data structures

ci_frac = zeros(2,K,K+1);
ci_yvalues_frac = zeros(2,K,K+1);
ci_check_frac = zeros(2,K,K+1);

%% Estimation

tol = 1e-3;

[params_simple_frac, optlikelihood_simple_frac] = estimation_simple_frac(data_frac,alpha,Ntot,T,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac);
[x0_frac, scale_frac] = initialguess_frac(params_simple_frac,lb_frac,ub_frac);
[params_mle_frac, error_mle_frac, scale_mle_frac, scale_error_frac, optlikelihood_frac, feasible_frac] ...
 = estimation_mle_frac_error(data_frac,C_frac,Ntot,alpha,T,x0_frac,scale_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);

if ~isempty(x0_frac_def)
    for n = 1:size(x0_frac_def,3)
        scale_frac_def = 10.^floor(log10(abs(x0_frac_def(:,:,n))));
        scale_frac_def(find(scale_frac_def == 0)) = 1;
        [params_mle_temp, error_mle_temp, scale_temp, scale_error_temp, optvalue_temp, feasible_temp]  ...
            = estimation_mle_frac_error(data_frac,C_frac,Ntot,alpha,T,x0_frac_def(:,:,n),scale_frac_def,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);
        
        if optvalue_temp < optlikelihood_frac
            params_mle_frac = params_mle_temp;
            error_mle_frac = error_mle_temp;
            scale_mle_frac = scale_temp;
            scale_error_frac = scale_error_temp;
            optlikelihood_frac = optvalue_temp;
            feasible_frac = feasible_temp;
        end              
    end
end

if nopt_simple_mle > 0 || nopt_random_mle > 0
    for n = 1:nopt_simple_mle
        [x0_frac, scale_frac] = initialguess_frac(params_simple_frac,lb_frac,ub_frac);
        [params_mle_temp, error_mle_temp, scale_temp, scale_error_temp, optvalue_temp, feasible_temp]  ...
            = estimation_mle_frac_error(data_frac,C_frac,Ntot,alpha,T,x0_frac,scale_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);

        if optvalue_temp < optlikelihood_frac
            params_mle_frac = params_mle_temp;
            error_mle_frac = error_mle_temp;
            scale_mle_frac = scale_temp;
            scale_error_frac = scale_error_temp;
            optlikelihood_frac = optvalue_temp;
            feasible_frac = feasible_temp;
        end     
    end
    for n = 1:nopt_random_mle
        [x0_frac, scale_frac] = initialguess_frac([],lb_frac,ub_frac);
        [params_mle_temp, error_mle_temp, scale_temp, scale_error_temp, optvalue_temp, feasible_temp]  ...
            = estimation_mle_frac_error(data_frac,C_frac,Ntot,alpha,T,x0_frac,scale_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,tol);

        if optvalue_temp < optlikelihood_frac
            params_mle_frac = params_mle_temp;
            error_mle_frac = error_mle_temp;
            scale_mle_frac = scale_temp;
            scale_error_frac = scale_error_temp;
            optlikelihood_frac = optvalue_temp;
            feasible_frac = feasible_temp;
        end     
    end
end

%% Confidence intervals
quant = chi2inv(1-alpha_q,1);
for i=1:K
    for j=1:K+1
        if ci_option(i,j) == 1

            ci_frac(1,i,j) = Inf;
            ci_frac(2,i,j) = -Inf;

            counter_1 = 0;
            counter_2 = 0;

            for n=1:nopt_ci+1
                [ci_frac_left_temp, ci_yvalues_frac_left_temp, ci_check_frac_left_temp] ...
                    = ci_left_endpoint_frac_error(data_frac,C_frac,Ntot,alpha,T,params_mle_frac,scale_mle_frac,optlikelihood_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,i,j,quant,tol);
                
                if ci_check_frac_left_temp == 0 || ci_frac_left_temp >= params_mle_frac(i,j)
                    end_cond = 0;  
                    while end_cond == 0
                        if mod(counter_1,2) == 0
                            [x0_frac, scale_frac] = initialguess_frac(params_mle_frac,lb_frac,ub_frac);
                        else
                            [x0_frac, scale_frac] = initialguess_frac([],lb_frac,ub_frac);
                        end
                        [ci_frac_left_temp, ci_yvalues_frac_left_temp, ci_check_frac_left_temp]...
                            = ci_left_endpoint_frac_error(data_frac,C_frac,Ntot,alpha,T,params_mle_frac,scale_mle_frac,optlikelihood_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,i,j,quant,tol);
    
                        if ci_check_frac_left_temp == 1 && ci_frac_left_temp <= params_mle_frac(i,j)
                            end_cond = end_cond+1;
                        end
                        counter_1 = counter_1+1;
                    end
                end
                
                if ci_frac_left_temp <= ci_frac(1,i,j)
                    ci_frac(1,i,j) = ci_frac_left_temp;
                    ci_yvalues_frac(1,i,j) = ci_yvalues_frac_left_temp;
                    ci_check_frac(1,i,j) = ci_check_frac_left_temp;
                end
    
                [ci_frac_right_temp, ci_yvalues_frac_right_temp, ci_check_frac_right_temp] ...
                    = ci_right_endpoint_frac_error(data_frac,C_frac,Ntot,alpha,T,params_mle_frac,scale_mle_frac,optlikelihood_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,i,j,quant,tol);
                
                if ci_check_frac_right_temp == 0 || ci_frac_right_temp <= params_mle_frac(i,j)
                    end_cond = 0;  
                    while end_cond == 0
                        if mod(counter_2,2) == 0
                            [x0_frac, scale_frac] = initialguess_frac(params_mle_frac,lb_frac,ub_frac);
                        else
                            [x0_frac, scale_frac] = initialguess_frac([],lb_frac,ub_frac);
                        end
                         [ci_frac_right_temp, ci_yvalues_frac_right_temp, ci_check_frac_right_temp] ...
                            = ci_right_endpoint_frac_error(data_frac,C_frac,Ntot,alpha,T,params_mle_frac,scale_mle_frac,optlikelihood_frac,Aineq_frac,bineq_frac,Aeq_frac,beq_frac,lb_frac,ub_frac,i,j,quant,tol);
    
                        if ci_check_frac_right_temp && ci_frac_right_temp >= params_mle_frac(i,j)
                            end_cond = end_cond+1;
                        end
                        counter_2 = counter_2+1;
                    end
                end
                
                if ci_frac_right_temp >= ci_frac(2,i,j)
                    ci_frac(2,i,j) = ci_frac_right_temp;
                    ci_yvalues_frac(2,i,j) = ci_yvalues_frac_right_temp;
                    ci_check_frac(2,i,j) = ci_check_frac_right_temp;
                end
            end
        end
    end
end

%% Data visualization
quant_norm = norminv(1-alpha_q/2);
B = eye(K,K-1);
if data_vis_option == 1
    [lambda_rel,d,lambda_1,A_rel] = params_extract_frac(params_mle_frac);
    for i=1:I
        for j=1:size(C_frac{i},2)
            figure;
            hold on
            c_j = B*C_frac{i}(:,j);
            n = size(T,2)*10;
            x = [0:T(end)/n:T(end)];
            f_mean = zeros(n,1);
            f_lb = zeros(n,1);
            f_ub = zeros(n,1);
            for l=1:n+1
                f_mean(l) = N(i,:)*expm(x(l)*A_rel)*c_j/sum(N(i,:)*expm(x(l)*A_rel));
                f_lb(l) = f_mean(l)-quant_norm*sqrt(sum(N(i,:))^(-1)*c_j'*covmatrix_frac(alpha(i,:),A_rel,lambda_rel,d,lambda_1,x(l))*c_j+c_j'*error_mle_frac^2*c_j);
                f_ub(l) = f_mean(l)+quant_norm*sqrt(sum(N(i,:))^(-1)*c_j'*covmatrix_frac(alpha(i,:),A_rel,lambda_rel,d,lambda_1,x(l))*c_j+c_j'*error_mle_frac^2*c_j);
            end
            plot(x,f_lb);
            plot(x,f_ub);
            patch([x fliplr(x)], [f_lb' fliplr(f_ub')], [0.95 0.95 0.95]);
            plot(x,f_mean,'Color',[0 0 0]);
            y = zeros(1,L);
            for r=1:R
                for l=1:L
                    y(l) = data_frac(i,:,l,r)*c_j;
                end
                scatter([0, T], [alpha(i,:)*c_j, y], 200, [0 0 0],'x','LineWidth',1);
            end
        end
    end
end
