function[grain_yield, straw_yield] = estimate_yield(Datam, RB, OPT_only_positive_yield)
    if nargin == 2
        OPT_only_positive_yield = false;
    end
    
    years_daily = Datam(1:24:end,1); %! Year for each day of the simulation (Datam is in h, hence only every 24th)
    years = unique(years_daily);
    years = years(2:end); % Ignore first year, which has no harvest
    yield_yearly_pred = zeros(length(years), size(RB, 2), 7);
    for i=1:length(years)
        yr = years(i);
        RB_curr_yr = RB(years_daily==yr, :, :);
        if OPT_only_positive_yield
            RB_curr_yr(RB_curr_yr < 0) = 0;
        end
        yield_yearly_pred(i, :, :) = sum(RB_curr_yr, 1);
    end
    grain_yield = sum(yield_yearly_pred(:, :, [4, 5]), 3); %! Grains are fruit (5) and C reserve (4) compartments
    straw_yield = sum(yield_yearly_pred(:, :, [1, 2, 7]), 3); %! Straw is leaves (1), sapwood (2), and dead leaves (7) compartments

end