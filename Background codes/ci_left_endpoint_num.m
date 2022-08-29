function [ci_lb,ci_yvalues_lb,ci_lb_check] = ci_left_endpoint_num(data,dead,C,N,T,mle,scale,optlikelihood,Aineq,bineq,Aeq,beq,lb,ub,i,j,quant,tol)

ci_lb_check = 0;

[Aineq,bineq,Aeq,beq,lb,ub] = addscaleconstraints_num(Aineq,bineq,Aeq,beq,lb,ub,scale);
mle = mle./scale;

options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-10,'TolCon',1e-10,'TolX',1e-10,'MaxFunEvals',50000,'MaxIter',10000,'ScaleProblem',true);
sol = fmincon(@(param) param(i,j),mle,Aineq,bineq,Aeq,beq,lb,ub,@(param) constraints_num(param,scale,C,N,data,dead,T,optlikelihood,quant),options);

ci_yvalues_lb = neg_log_likelihood_num(sol,scale,C,N,data,dead,T);
sol = sol.*scale;
ci_lb = sol(i,j);

if (abs(ci_yvalues_lb - optlikelihood - quant)/(optlikelihood+quant)<tol) || (ci_lb == lb(i,j) && ci_yvalues_lb < optlikelihood + quant + tol)
    ci_lb_check = 1;
end

end

function [c,ceq] = constraints_num(param,scale,C,N,data,dead,T,optvalue,quant)
    c = neg_log_likelihood_num(param,scale,C,N,data,dead,T) - optvalue - quant;
    ceq = 0;
end