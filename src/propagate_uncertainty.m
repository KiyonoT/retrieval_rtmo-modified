function uncertainty = propagate_uncertainty(paramsout, measurement, tab, angles, irr_meas, fixed, sensor)    

    tab.value = paramsout;
    iparams = tab.include;
    Sa = tab.uncertainty(iparams);
    
    %% here modification of parameters should occur
    tab = helpers.modify_tab_parameters(tab);
        
    %% save prior (x0) values
    tab.x0 = tab.value;
    
    tab.include = true(numel(tab.include), 1);
    J = helpers.clc_jacobian(-999, measurement, tab, angles, irr_meas, fixed, sensor);
    J = J(measurement.i_fit,:);

    % propagation of uncertainty in refl to uncertainty in variables?
    meas_std_fit = measurement.std(measurement.i_fit);
    meas_std_fit(isinf(meas_std_fit)) = 0;
%     std = abs(pinv(J)*meas_std_fit); % Moore-Penlose solution
    K = J(:,iparams);
    G = diag(Sa.^2)*K.' / (K*diag(Sa.^2)*K.' + diag(meas_std_fit.^2)); % Gain matrix (Rodgers 2000 [Eq 2.45])
    std = NaN(1, length(tab.uncertainty));
    std(iparams) = abs(G*meas_std_fit); % MAP solution
    
    % posterior parameter uncertainty
    K = K(meas_std_fit~=0, :);
    meas_std_fit = meas_std_fit(meas_std_fit~=0);
    post = NaN(1, length(tab.uncertainty));
    post(iparams) = diag((diag(Sa.^-2) + K.'*diag(meas_std_fit.^-2)*K).^-0.5);

    % propagation of uncertainty in variables to uncertainty in reflectance
    in_covar_mat = diag(tab.uncertainty);  % we do not have and do not need covariances
    out_covar_mat = J * in_covar_mat * J';
    std_refl = diag(out_covar_mat);  % on diagonal var are located

    uncertainty.J = J;
    uncertainty.gained_noise = std;
    uncertainty.std_params = post;
    uncertainty.std_refl = std_refl;
    
%     figure(c*1000)
% %     plot(measurement.wl(measurement.i_fit), std_j)
%     errorbar(measurement.wl(measurement.i_fit), reflSAIL(measurement.i_fit), std_j)
%     ylim([0, 1])
end