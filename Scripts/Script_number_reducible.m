%% Background information

% This is a script for cell number data without measurement error, for a
% branching process model with reducible switching dynamics. See Appendix B
% of the accompanying paper for a discussion of such a model.

% The script implements the statistical model (5) of Section 4.2 of the
% accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error
% term. The data is artifical and generated from the model structure shown
% in Figure 10 of the accompanying paper. The parameters are b_1 = 0.6, b_2
% = 1.0, b_3 = 0.8, d_1 = 0.3, d_2 = 0.5, d_3 = 0.6, nu_(12) = 0.02,
% nu_(21) = 0.01 and nu_(23) = 0.003. All other switching rates are 0.

% In the script, the parameters of the branching process model from Section
% 2 of the accompanying paper are arranged in a K x (K+1) matrix as
% follows:

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

% params_simple_num is a K x (K+1) matrix which holds the estimates for the
% model parameters under the deterministic population model presented in
% Appendix C of the accompanying paper. Arrangement of parameters is as
% described above.
 
% params_mle_num is a K x (K+1) matrix which holds the estimates for
% the model parameters under the statistical model (5) of Section 4.2, with
% no measurement error term.

% If confidence intervals are requested, ci_num is a 2 x K x (K+1)
% tensor, where for each i = 1,...,K and j = 1,...,K+1, ci_frac(:,i,j)
% holds the lower and upper interval bounds for parameter (i,j), with the
% arrangement of parameters described at the top of the script.

% If requested, the script also produces a visualization of how well the
% model fits the data. See further explanation below.

%% Input parameters

% A random seed can be fixed for reproducibility.
rndseed = 1; rng(rndseed);

% K is the number of phenotypes in the model.
K = 3;
% I is the number of initial conditions used in the cell line experiments.
I = 3;
% R is the number of experimental replicates performed.
R = 3;
% T is the row vector of experimental timepoints.
T = [1:1:6];
% L is the number of experimental timepoints.
L = size(T,2);

% data_num is an I x K x L x R tensor that holds the experimental data.
% Entry number (i,k,l,r) is the number of cells of type-k in the r-th
% replicate of the experiment started by the i-th initial condition and
% ended at the l-th timepoint.
% Here, we have manually entered artifical data that was simulated from the
% branching process model of Section 2 of the accompanying paper. The
% parameters of the model are as stated at the top of the script.
data_num = zeros(I,K,L,R);
data_num(1,1,:,1) = [13206,17468,22838,30621,41498,54904];
data_num(1,2,:,1) = [282,845,1921,3886,6874,14029];
data_num(1,3,:,1) = [0,2,5,7,21,79];
data_num(2,1,:,1) = [140,459,926,1862,3596,6382];
data_num(2,2,:,1) = [16205,26481,42224,70730,115321,185091];
data_num(2,3,:,1) = [44,110,222,452,857,1647];
data_num(3,1,:,1) = [0,0,0,0,0,0];
data_num(3,2,:,1) = [0,0,0,0,0,0];
data_num(3,3,:,1) = [12411,14992,17652,21730,26934,32931];
data_num(1,1,:,2) = [12964,17762,23492,30633,41268,54260];
data_num(1,2,:,2) = [324,833,2020,3949,6809,13239];
data_num(1,3,:,2) = [0,3,7,13,52,64];
data_num(2,1,:,2) = [171,453,925,1933,3523,6331];
data_num(2,2,:,2) = [16301,26620,43428,68892,115914,189285];
data_num(2,3,:,2) = [46,131,201,413,874,1709];
data_num(3,1,:,2) = [0,0,0,0,0,0];
data_num(3,2,:,2) = [0,0,0,0,0,0];
data_num(3,3,:,2) = [11918,14671,18975,23094,27270,33453];
data_num(1,1,:,3) = [13291,17686,23613,30331,40761,54108];
data_num(1,2,:,3) = [263,834,1895,3995,7596,12313];
data_num(1,3,:,3) = [0,1,9,22,45,109];
data_num(2,1,:,3) = [138,415,965,1906,3657,6657];
data_num(2,2,:,3) = [16492,25582,44012,67864,117040,188723];
data_num(2,3,:,3) = [44,115,270,444,1100,1709];
data_num(3,1,:,3) = [0,0,0,0,0,0];
data_num(3,2,:,3) = [0,0,0,0,0,0];
data_num(3,3,:,3) = [12230,14764,18504,21139,27339,33157];

% The user has the option of supplying an I x 1 x L x R tensor with the
% number of dead cells at each time point. In that case, the estimation is
% performed on the agumented model shown in Figure 6 of the accompanying
% paper, where a new state representing dead cells is added.
% Note that the augmented model does not take into account clearance of
% dead cells. See Section 7.5 of the accompanying paper.
dead = [];

% N is an I x K matrix where the i-th row contains the starting number of
% cells of each type under the i-th initial condition.
N = [10^4,0,0;0,10^4,0;0,0,10^4];

% C_num is a cell array, where C_num{i} is a K x J_i matrix which
% allows the user to reduce the experimental data under the i-th initial
% condition. In this case, C_num{3} is set to be the 3rd unit column
% vector, to account for the fact that type-3 cells do not switch over to
% type-1 or type-2. See Appendix B of the accompanying paper.
C_num = cell(I);
for i=1:I-1
    C_num{i} = eye(K,K);
end
C_num{I} = [0;0;1];

% lb_num is a K x (K+1) matrix which holds lower bounds for the branching
% process model parameters. The arrangement of parameters is as described
% at the top of the script. For example, the first column holds lower
% bounds for the birth rates of the phenotypes.
lb_num = [0,-4,0,0;0,-4,0,0;0,-4,0,0];
% ub_num is a K x (K+1) matrix which holds upper bounds for the branching
% process model parameters.
ub_num = [2,4,1,1;2,4,1,1;2,4,1,1];
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
Aeq_num = zeros(K,K+1,3);
Aeq_num(:,:,1) = [0,0,0,1;0,0,0,0;0,0,0,0];
Aeq_num(:,:,2) = [0,0,0,0;0,0,0,0;0,0,1,0];
Aeq_num(:,:,3) = [0,0,0,0;0,0,0,0;0,0,0,1];
% beq_num is an n x 1 vector, where n is the number of linear equality
% constraints. For each m = 1,...,n, beq_num(m) holds the coefficient of
% the right-hand side of the m-th equality constraint.
beq_num = [0;0;0];
% In this script, the switching rates nu_(13), nu_(31) and nu_(32) are
% constrained to be zero. To impose these constraints, the coefficients of
% nu_(13), nu_(31) and nu_(32) in Aeq_num are set to 1 respectively, and % 
% beq_num is set to a three-dimensional column vector of zeros.

% x0_num_def is a K x (K+1) x n tensor, which can be used to supply
% initial guesses for the MLE optimization problem. Here, n is the number
% of distinct initial guesses the user wishes to supply. The arrangement of
% parameters is as described at the top of the script.
x0_num_def = [];

% By default, the MLE optimization problem is solved once, using an initial
% guess based on simple deterministic parameter estimates. Here, the user
% is given the option to request that the optimization problem is solved
% for multiple initial guesses, each based on the simple deterministic
% estimate. See Appendix C of the accompanying paper.
nopt_simple_mle = 0;
% Here, the user is given the option to request that the MLE optimization
% problem for is solved for multiple initial guesses, where the initial
% guesses are chosen in a random fashion as described in Appendix D of the
% accompanying paper.
nopt_random_mle = 0;

% ci_option is a K x (K+1) matrix holding zeros and ones. Setting
% ci_option(i,j) = 1 will compute a confidence interval for parameter (i,j)
% of the branching process model, where the arrangement of parameters is
% as described at the top of the script.
% In this script, a confidence interval is requested for all parameters
% except the switching rates known to be zero.
ci_option = ones(K,K+1);
ci_option(1,4) = 0;
ci_option(3,3) = 0;
ci_option(3,4) = 0;
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

ci_num = zeros(2,K,K+1);
ci_yvalues_num = zeros(2,K,K+1);
ci_check_num = zeros(2,K,K+1);

%% Estimation

tol = 1e-3;

params_simple_num = estimation_simple_num(data_num,N,T);   
[x0_num, scale_num] = initialguess_num(dead,N,T,params_simple_num,lb_num,ub_num);
[params_mle_num, scale_mle_num, optlikelihood_num, feasible_num] ...
 = estimation_mle_num(data_num,dead,C_num,N,T,x0_num,scale_num,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,tol);

if ~isempty(x0_num_def)
    scale_num_def = 10.^floor(log10(abs(x0_num_def(:,:,n))));
    scale_num_def(find(scale_num_def == 0)) = 1;
    for n = 1:size(x0_num_def,3)
        [params_mle_temp, scale_temp, optvalue_temp, feasible_temp] ...
            = estimation_mle_num(data_num,dead,C_num,N,T,x0_num_def(:,:,n),scale_num_def,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,tol);
        
        if optvalue_temp < optlikelihood_num
            params_mle_num = params_mle_temp;
            scale_mle_num = scale_temp;
            optlikelihood_num = optvalue_temp;
            feasible_num = feasible_temp;
        end              
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

%% Confidence intervals
quant = chi2inv(1-alpha_q,1);

for i=1:K
    for j=1:K+1
        if ci_option(i,j) == 1

            ci_num(1,i,j) = Inf;
            ci_num(2,i,j) = -Inf;

            counter_1 = 0;
            counter_2 = 0;
            
            for n=1:nopt_ci+1
                end_cond = 0;   
                while end_cond == 0
                    if mod(counter_1,2) == 0
                        [x0_num, scale_num] = initialguess_num(dead,N,T,params_mle_num,lb_num,ub_num);
                    else
                        [x0_num, scale_num] = initialguess_num(dead,N,T,[],lb_num,ub_num);
                    end
                    [ci_num_left_temp, ci_yvalues_num_left_temp, ci_check_num_left_temp] ...
                        = ci_left_endpoint_num(data_num,dead,C_num,N,T,params_mle_num,scale_mle_num,optlikelihood_num,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,i,j,quant,tol);

                    if ci_check_num_left_temp == 1 && ci_num_left_temp <= params_mle_num(i,j)
                        end_cond = end_cond+1;
                    end
                    counter_1 = counter_1+1;
                end

                if ci_num_left_temp <= ci_num(1,i,j)
                    ci_num(1,i,j) = ci_num_left_temp;
                    ci_yvalues_num(1,i,j) = ci_yvalues_num_left_temp;
                    ci_check_num(1,i,j) = ci_check_num_left_temp;
                end

                end_cond = 0;   
                while end_cond == 0
                    if mod(counter_2,2) == 0
                        [x0_num, scale_num] = initialguess_num(dead,N,T,params_mle_num,lb_num,ub_num);
                    else
                        [x0_num, scale_num] = initialguess_num(dead,N,T,[],lb_num,ub_num);
                    end
                    [ci_num_right_temp, ci_yvalues_num_right_temp, ci_check_num_right_temp] ...
                        = ci_right_endpoint_num(data_num,dead,C_num,N,T,params_mle_num,scale_mle_num,optlikelihood_num,Aineq_num,bineq_num,Aeq_num,beq_num,lb_num,ub_num,i,j,quant,tol);

                    if ci_check_num_right_temp && ci_num_right_temp >= params_mle_num(i,j)
                        end_cond = end_cond+1;
                    end
                    counter_2 = counter_2+1;
                end

                if ci_num_right_temp >= ci_num(2,i,j)
                    ci_num(2,i,j) = ci_num_right_temp;
                    ci_yvalues_num(2,i,j) = ci_yvalues_num_right_temp;
                    ci_check_num(2,i,j) = ci_check_num_right_temp;
                end
            end
        end
    end
end

%% Data visualization
quant_norm = norminv(1-alpha_q/2);
if data_vis_option == 1
    [b,A] = params_extract_num(params_mle_num,[]);
    for i=1:I
        for j=1:size(C_num{i},2)
            figure;
            hold on
            c_j = C_num{i}(:,j);
            n = size(T,2)*10;
            x = [0:T(end)/n:T(end)];
            f_mean = zeros(n,1);
            f_lb = zeros(n,1);
            f_ub = zeros(n,1);
            for l=1:n+1
                f_mean(l) = N(i,:)*expm(x(l)*A)*c_j;
                f_lb(l) = f_mean(l)-quant_norm*sqrt(c_j'*covmatrix_num(N(i,:),A,b,x(l),[])*c_j);
                f_ub(l) = f_mean(l)+quant_norm*sqrt(c_j'*covmatrix_num(N(i,:),A,b,x(l),[])*c_j);
            end
            plot(x,f_lb);
            plot(x,f_ub);
            patch([x fliplr(x)], [f_lb' fliplr(f_ub')], [0.95 0.95 0.95]);
            plot(x,f_mean,'Color',[0 0 0]);
            y = zeros(1,L);
            for r=1:R
                for l=1:L
                    y(l) = data_num(i,:,l,r)*c_j;
                end
                scatter([0, T], [N(i,:)*c_j, y], 200, [0 0 0],'x','LineWidth',1);
            end
        end
    end
end