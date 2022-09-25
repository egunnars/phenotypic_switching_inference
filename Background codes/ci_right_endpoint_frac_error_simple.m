function [ci_ub,ci_yvalues_ub,ci_ub_check] = ci_right_endpoint_frac_error_simple(data,C,Ntot,alpha,T,mle,scale,optlikelihood,Aineq,bineq,Aeq,beq,lb,ub,i,j,quant,tol)

ci_ub_check = 0;

K = size(alpha,2);
[Aineq,bineq,Aeq,beq,lb,ub] = addscaleconstraints_frac(Aineq,bineq,Aeq,beq,lb,ub,scale);
mle = mle./scale;

Aineq_new = [Aineq,zeros(size(Aineq,1),1)];
if ~isempty(Aeq)
   Aeq_new = [Aeq,zeros(size(Aeq,1),1)];
else
   Aeq_new = Aeq;
end
scale_new = [reshape(scale,1,K*(K+1)),0.01];
lb_new = [reshape(lb,1,K*(K+1)),0];
ub_new = [reshape(ub,1,K*(K+1)),1/scale_new(K*(K+1)+1)];
mle_new = [reshape(mle,1,K*(K+1),1),0.01/scale_new(K*(K+1)+1)];

options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',50000,'MaxIter',10000,'ScaleProblem',true);
idx = sub2ind([K K+1],i,j);
sol = fmincon(@(param) -param(idx),mle_new,Aineq_new,bineq,Aeq_new,beq,lb_new,ub_new,@(param) constraints_frac(param,scale_new,C,Ntot,alpha,data,T,optlikelihood,quant),options);

ci_yvalues_ub = neg_log_likelihood_frac_error_simple(sol,scale_new,C,Ntot,alpha,data,T);
sol = reshape(sol(1:K*(K+1)),K,K+1);
sol = sol.*scale;
ci_ub = sol(i,j);

if (abs((ci_yvalues_ub - optlikelihood - quant)/(optlikelihood+quant))<tol) || (ci_ub == ub(i,j) && ci_yvalues_ub < optlikelihood + quant + tol)
    ci_ub_check = 1;
end

end

function [c,ceq] = constraints_frac(param,scale,C,Ntot,alpha,data,T,optvalue,quant)
    c = neg_log_likelihood_frac_error_simple(param,scale,C,Ntot,alpha,data,T) - optvalue - quant;
    ceq = 0;
end