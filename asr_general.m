function [estimates, convergence_iter, likelihood_per_iter] = asr_general(data, lambda, epsilon, max_iter)
    
    % data is an mx(n+1) matrix, where m is the total number of comparisons
    % and n is the number of items. Each row of this matrix is of the
    % following format: 
    %   data(i,1) = winner of the i^th comparison
    %   forall j in {2,..,n+1}, data(i,j) = 1 if element j-1 was a part of the i^th comparison
    
    % lambda is the regularization parameter
    % epsilon is the threshold value for convergence
    
    sets = data(:,2:end);
    wins = data(:,1);
    unique_sets = unique(sets,'rows');

    n = size(sets,2);
    d = zeros(1,n); % degrees of element vertices in the comparison graph
    for i = 1:size(sets,1)
       m = sum(sets(i,:)); 
       d(logical(sets(i,:)))=d(logical(sets(i,:))) + 1/m;
    end

    for i = 1:size(unique_sets,1)
       d(logical(unique_sets(i,:))) = d(logical(unique_sets(i,:)))+lambda;
    end
    
    P = zeros(n,n);
    for i = 1:size(unique_sets,1)
       current_set = unique_sets(i,:);
       m = sum(current_set);
       for j = 1:size(current_set,2)
          if current_set(j) == 0
              continue;
          end
          for k = 1:size(current_set,2)
             if current_set(k) == 0
                 continue;
             end
             P(j,k) = P(j,k) + lambda/m;
          end
       end
    end
    for i = 1:size(sets,1)
       current_set = sets(i,:);
       current_win = wins(i);
       m = sum(current_set);
       current_ele = [1:n];
       current_ele = current_ele(logical(current_set));
       for j = 1:size(current_ele,2)
           P(current_ele(j),current_win) = P(current_ele(j),current_win) + 1/m;
       end
    end
    for row = 1:n
        P(row,:) = P(row,:)./d(row);
    end

    estimates = ones(1,n)/n;
    estimates_old = estimates;
    likelihood_per_iter = zeros(max_iter,1);
    convergence_iter = n;

    for iter = 1:max_iter

        estimates = estimates*P;
        estimates = estimates/sum(estimates);
        
        likelihood_per_iter(iter) = loglikelihood(sets, wins, (estimates./d)/sum(estimates./d));

        if(norm(estimates-estimates_old,2)< epsilon*norm(estimates_old,2))
            convergence_iter = iter;
            break;
        end
        if(iter == max_iter)
            disp('Accelerated spectral ranking failed to converge\n');
        end

        estimates_old = estimates;
    end
    
    estimates = (estimates./d)/sum(estimates./d);
    likelihood_per_iter(convergence_iter:end) = likelihood_per_iter(convergence_iter);

end