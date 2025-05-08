function[B_IO] = load_fertilizer(B_IO, datenum_start, num_days, crop_data, manure_data_path, OPT_Use_Fertilizer)
    if nargin == 4
        % Assume that fertilizer to be used by default
        OPT_Use_Fertilizer = true;
    end

    Time = datetime( (datenum_start:datenum_start+num_days)', "convertfrom", "datenum");
    fert_dates = timetable(Time, zeros(size(Time)), zeros(size(Time)), zeros(size(Time)), zeros(size(Time)),...
        'VariableNames', {'Nrate', 'Prate', 'Krate', 'ManF'});
    
    if OPT_Use_Fertilizer
        % Select only entries where N is applied
        Nindices = ~isnat(crop_data.spring2_n_date);
        Ndates = crop_data.spring2_n_date(Nindices);
        Nrates = crop_data.spring2_n_kg_ha(Nindices) * 0.1; % 1 kg/ha = 0.1 g/m2
    
        for i=1:length(Ndates)
            date_curr = Ndates(i);
            % Make sure the date is in simulated time range
            if ~isempty(fert_dates(date_curr,:))
                fert_dates(date_curr,'Nrate') = {Nrates(i)};
            end
        end
    
        %%%% Manure data
    
        % Only select relevant columns
        opts = detectImportOptions(manure_data_path);
        opts.SelectedVariableNames = {'harvest_year', 'manure_type', 'C_kg_ha', 'C_N', 'P_kg_ha', 'K_kg_ha', 'DM_percent', 'C_percent_DM'};
        opts = setvartype(opts, {'C_kg_ha', 'C_N', 'P_kg_ha', 'K_kg_ha', 'DM_percent', 'C_percent_DM'}, 'double');
        manure_data = readtable(manure_data_path, opts);
    
        % Remove entries with invalid data
        manure_idx = ~isnan(manure_data.harvest_year)...
            & strcmp(manure_data.manure_type, 'FYM')...
            & ~isnan(manure_data.C_kg_ha)...
            & ~isnan(manure_data.P_kg_ha)...
            & ~isnan(manure_data.K_kg_ha)...
            & ~isnan(manure_data.C_N);
        manure_data = manure_data(manure_idx,:);
    
        % P and K fractions
        manure_data.C_P = manure_data.C_kg_ha ./ manure_data.P_kg_ha;
        manure_data.C_K = manure_data.C_kg_ha ./ manure_data.K_kg_ha;
    
        % For now take mean values across different years for manure composition,
        % even though it does vary substantially
        B_IO.N_Man = mean(manure_data.C_N);
        B_IO.P_Man = mean(manure_data.C_P);
        B_IO.K_Man = mean(manure_data.C_K);
        B_IO.Lig_fr_Man = 0.1; % Ligning fraction for cow manure is typically between 6 and 14%
    
        % Mean C content
        manure_C_content = mean(manure_data.DM_percent .* manure_data.C_percent_DM) * 1e-4;  % 1e-4 because values are in per cent
    
        % Select rows where manure is applied
        crop_data_manure = crop_data(strcmp(crop_data.fym_factor_level, "FYM"), :);
        Mdates = datetime(crop_data_manure.fym_date);
        Mrates = crop_data_manure.fym_applied * manure_C_content * 100; % 1 t/ha = 100 g/m2
    
        for i=1:length(Mdates)
            date_curr = Mdates(i);
            % Make sure the date is in simulated time range
            if ~isempty(fert_dates(date_curr,:))
                fert_dates(date_curr,"ManF") = {Mrates(i)};
            end
        end
    end
    
    B_IO.FertN = fert_dates.Nrate;
    B_IO.FertP = fert_dates.Prate;
    B_IO.FertK = fert_dates.Krate;
    B_IO.ManF = fert_dates.ManF;