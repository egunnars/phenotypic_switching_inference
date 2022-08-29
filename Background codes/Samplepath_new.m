function [num, num_dead] = Samplepath_new(params_num,N,T)

K = size(N,2);
num = N;

params = params_num;
params(:,2) = params(:,1)-params(:,2);

a = sum(params,2)';
norm_params = params./a';
sum_norm_params = zeros(K,K+1);
for i=1:K
    for j=1:K+1
        sum_norm_params(i,j) = sum(norm_params(i,1:j));
    end
end

t = 0;
num_dead = 0;

while sum(num) > 0 && t < T
    totalrates = num.*a;
    sumtotalrates = sum(totalrates);

    cumsum = 0;
    t = t - 1/(sum(totalrates))*log(rand);
    
    u = rand;
    for i=1:K
        cumsum = cumsum + totalrates(i);
        if u < cumsum/sumtotalrates
            numtype = i;
            break;
        end
    end

    u = rand;
    for j=1:K+1
        if u < sum_norm_params(numtype,j)
           numevent = j;
           break;
        end
    end
    
    if numevent == 1
        num(numtype) = num(numtype)+1;
        else if numevent == 2
            num(numtype) = num(numtype)-1;
            num_dead = num_dead+1;
            else if numevent-2 < numtype
                num(numtype) = num(numtype)-1;
                num(numevent-2) = num(numevent-2)+1;
            else
                num(numtype) = num(numtype)-1;
                num(numevent-1) = num(numevent-1)+1;
            end
        end
    end    
end