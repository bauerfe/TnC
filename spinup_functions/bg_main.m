function bg_main(initial_spinup_results_path, initial_spinup_id, relevant_strip, treatment_descriptor, path_deposition_data, ...
    manure_data_path, Nyears_epoch, N_epochs_max, rtol, atol, OPT_Use_Fertilizer, OPT_Plot)

    % Set up output directory
    bg_descriptor = [sprintf('%d', relevant_strip) '_' treatment_descriptor];
    if length(bg_descriptor) > 0
        run_identifier = [initial_spinup_id '_' bg_descriptor];
    else
        run_identifier = initial_spinup_id;
    end
    bg_save_folder = ['results' filesep 'bg_spinup' filesep run_identifier];
    mkdir(bg_save_folder);

    % Load input data from initial spinup
    load(initial_spinup_results_path);
    
    % Preprocessing - Function arguments are from initial spinup
    [Se, Se_fc, Psi_s, Ts, V, VT] = bg_preprocess(Ta, Tdp, O, V, Soil_Param, Phy, SPAR, Bio_Zs);

    % The actual spinup
    [B, R_litter, R_litter_sur, R_microbe, R_bacteria, R_ew, VOL, BfixN, Min_N, Min_P, RmycAM, RmycEM, ...
        N2flx, NH4_Uptake, NO3_Uptake, P_Uptake, K_Uptake, LEAK_NH4, LEAK_NO3, LEAK_P, LEAK_K, LEAK_DOC, ...
        LEAK_DON, LEAK_DOP, Lk] = bg_spinup(Lat, Lon, N_epochs_max, Nyears_epoch, ISOIL_H+ISOIL_L, Zbio, rsd, ...
        Pcla, Psan, PHs, Zs, Ts, Ta, Psi_s, Se, Se_fc, V, VT, T_L, T_H, Lk, RexmyI, Ccrown, B_H+B_L, rtol, atol, bg_save_folder, ...
        dateNum1, crop_data, manure_data_path, OPT_Use_Fertilizer)
    
    % Plot results of final epoch (states should not change significantly anymore)
    if OPT_Plot
        bg_plot(B, R_litter, R_litter_sur, R_microbe, R_bacteria, R_ew, VOL, BfixN, Min_N, Min_P, RmycAM, RmycEM, ...
            N2flx, NH4_Uptake, NO3_Uptake, P_Uptake, K_Uptake, LEAK_NH4, LEAK_NO3, LEAK_P, LEAK_K, LEAK_DOC, ...
            LEAK_DON, LEAK_DOP, Lk, Zbio, rsd);
    end

end