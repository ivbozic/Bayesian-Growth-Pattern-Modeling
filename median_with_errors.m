function [median, low, high] = median_with_errors(a, low_percentile, high_percentile)

a_sorted = sort(a);
n = length(a_sorted);
if mod(n, 2) == 1
    median = a_sorted((n + 1) / 2);
else
    median = (a_sorted(n / 2) + a_sorted(n / 2 + 1)) / 2;
end
low = a_sorted(ceil(low_percentile * (n + 1)));
high = a_sorted(floor(high_percentile * (n + 1)));

    
end

