function [path_curr reached_equilibrium] = coupled_model_main(relevant_strip, descriptor_full, ...
    root_dir, output_dir, N_max_epochs, rtol, atol, compare_vars, prev_data_path, prev_bg_path, ...
    prova_file_path, initial_parameters_path, bg_parameters_path)

    prev_data = load(prev_data_path);
    if prev_bg_path=="none"
        prev_bg = prev_data.P(end,:);
    else
        prev_bg = load(prev_bg_path).B_end;
    end
    
    reached_equilibrium = false;
    
    for epoch=1:N_max_epochs
        disp(sprintf('Epoch %02d of %d', epoch, N_max_epochs))
        
        % Actual spinup
        [path_curr] = coupled_model_spinup(relevant_strip, descriptor_full, prev_data, prev_bg, ...
            root_dir, output_dir, prova_file_path, initial_parameters_path, bg_parameters_path, epoch)

        new_data = load(path_curr);

        % Test whether equilibrium has been reached
        reached_equilibrium = is_difference_inside_tol(prev_data, new_data, rtol, atol, compare_vars);
        if reached_equilibrium
            disp('    Reached dynamic equilibrium.')
            break;
        else
            prev_data = new_data;            
            prev_bg = new_data.P(end,:);
        end
    end

    if ~reached_equilibrium
        disp(sprintf('After %d epochs, dynamic equilibrium was not reached', epoch))
    end
end

function reached_equilibrium = is_difference_inside_tol(prev_data, new_data, rtol, atol, compare_vars)
    reached_equilibrium = true;
    
    for varname=compare_vars
        prev_val = getfield(prev_data, varname);
        new_val = getfield(new_data, varname);
        diff = abs(new_val - prev_val);
        % if any( (diff > atol) & (diff > rtol * new_val) )
        % For rtol take mean abs. diff, for atol take max value,
        % i.e. mean absolute error hass to be below rtol, and max error below atol
        % Still use `any` for 3D-arrays
        eps = 1e-10;
        prev_val(abs(prev_val) < eps) = eps;
        mean_rel_diff = mean(diff / prev_val);
        if any(diff > atol) | any(mean_rel_diff > rtol)
            reached_equilibrium = false;
            break
        end
    end
end  