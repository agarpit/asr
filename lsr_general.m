function [estimates, convergence_iter, likelihood_per_iter] = lsr_general(data, lambda, epsilon, max_iter)

    sets = data(:,2:end);
    wins = data(:,1);
    unique_sets = unique(sets,'rows');

    n = size(sets,2);
    P = zeros(n, n);

    for d = 1:size(unique_sets,1)
       current_set = unique_sets(d,:);
       m = sum(current_set);
       for i = 1:size(current_set,2)
          if current_set(i) == 0
              continue;
          end
          for j = 1:size(current_set,2)
             if current_set(j) == 0 || j == i
                 continue;
             end
             P(j,i) = P(j,i) + lambda/m;
          end
       end
    end
    
    for d = 1:size(data,1)
        current_set = sets(d,:);
        current_winner = wins(d);
        m = sum(current_set);
        for i = 1:size(current_set,2)
            if current_set(i) == 0 || i == current_winner
                continue;
            end
            P(i,current_winner) = P(i,current_winner) + 1/m;
        end
    end
    
    eps = max(sum(P,2));
    P = P - diag(sum(P,2));
    P = eye(n,n) + (1/eps)*P;

    estimates = ones(1,n)/n;
    estimates_old = estimates;
    likelihood_per_iter = zeros(max_iter,1);
    convergence_iter = n;

    for iter = 1:max_iter

        estimates = estimates*P;
        estimates = estimates/sum(estimates);

        likelihood_per_iter(iter) = loglikelihood(sets,wins,estimates);

        if(norm(estimates-estimates_old,2)<epsilon*norm(estimates_old,2))
            convergence_iter = iter;
            break;
        end
        if(iter==max_iter)
            disp('Luce spectral ranking failed to converge\n');
        end

        estimates_old = estimates;
    end
    
    likelihood_per_iter(convergence_iter:end) = likelihood_per_iter(convergence_iter);
    
end