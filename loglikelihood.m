function [ll] = loglikelihood(sets, wins, estimates)


num_sets = size(sets, 1);
n = size(sets,2);
ll = 0;

for i = 1:num_sets
    current_set = sets(i,:);
    current_win = wins(i);
    current_ele = [1:n];
    current_ele = current_ele(logical(current_set));
    
    x = 0;
    y = 0;
    if estimates(current_win) == 0
        x = 0;
    else
        x = log(estimates(current_win));
    end
    if sum(estimates(current_ele)) == 0
        y = 0;
    else
        y = log(sum(estimates(current_ele)));
    end
    ll = ll + (x-y);

end