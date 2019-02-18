function sample_value = sample_continuous(grid_points, grid_values)

cum_int = cumtrapz(grid_points, grid_values);
n_points = length(cum_int);

rand_val = rand(1) * cum_int(n_points);

left = 1;
right = n_points;
left_val = cum_int(left);
right_val = cum_int(right);

while right - left > 1
    mid = round((left + right) / 2);
    mid_val = cum_int(mid);
    if rand_val <= mid_val
        right = mid;
        right_val = mid_val;
    else
        left = mid;
        left_val = mid_val;
    end
end

delta_x = grid_points(right) - grid_points(left);
left_y = grid_values(left);
right_y = grid_values(right);
delta_y = right_y - left_y;
slope = delta_y / delta_x;
if abs(slope) <= 1e-10
    sample_value = grid_points(left) + delta_x / 2;
else
    sample_value = grid_points(left) + (sqrt(left_y^2 + 2 * slope * (rand_val - left_val)) - left_y) / slope;
end

end
