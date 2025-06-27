function [path_curr reached_equilibrium] = coupled_model_main(relevant_strip, descriptor_full, ...
    root_dir, output_dir, N_max_epochs, atol, rtol, compare_vars, prev_data_path, prev_bg_path, ...
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
        reached_equilibrium = is_difference_inside_tol(prev_data, new_data, atol, rtol, compare_vars);
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

function reached_equilibrium = is_difference_inside_tol(prev_data, new_data, atol, rtol, compare_vars)
    reached_equilibrium = true;
    
    for varname=compare_vars
        prev_val = getfield(prev_data, varname);
        new_val = getfield(new_data, varname);
        % % Where `prev_val` is greater than `atol`, take a relative difference
        % % and compare to `rtol`. Otherwise compare absolute difference to `atol`
        % small_vals = abs(prev_val) < atol;
        % abs_diff = abs(new_val - prev_val);
        % relative_diff = abs_diff(~small_vals) ./ abs(prev_val(~small_vals));
        % if any(relative_diff > rtol) | any(abs_diff(small_vals) > atol)
        rmse = sqrt( mean ((prev_val - new_val).^2) );
        if any(rmse > rtol * std(prev_val) + atol)  % Use `any` for arrays with >= 3 dim
            reached_equilibrium = false;
            break
        end
    end
end  