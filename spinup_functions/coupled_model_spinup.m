function save_path = coupled_model_spinup(relevant_strip, descriptor_full, prev_data, prev_bg, root_dir, ...
    output_dir, prova_file_path, initial_parameters_path, bg_parameters_path, epoch)
    %! `descriptor_full` is just passed to be saved by `save` function
    %! `relevant_strip`, prev_data, prev_bg will be used in PARAM_IC
    
    if nargin == 9
        epoch=1;
    elseif nargin ~= 10
        error('Function `coupled_model_spinup` got %d arguments. Should have received 7 or 8', nargin);        
    end
    
    % Run truncated prova file to perform standard initialization
    run(prova_file_path);
    OPT_BG = 1;
    OPT_Use_Fertilizer = 1;
    PARAM_IC = bg_parameters_path;
    % Run the simulation
    MAIN_FRAME

    mkdir(output_dir);
    save_path = determine_save_path(output_dir, epoch);
    save(save_path, '-regexp', '^(?!(prev_bg|prev_data)$).'); % Save data, but without saving previous data again
    disp(['Saved data under ' save_path])
end

function save_path = determine_save_path(output_dir, epoch)
    while true
        save_path_curr = [output_dir filesep sprintf('epoch_%02d', epoch) '.mat'];
        if isfile(save_path_curr) % File already exists
            epoch = epoch + 1;
        else
            save_path = save_path_curr;
            break;
        end
    end
end
    