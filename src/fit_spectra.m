function results = fit_spectra(measurement, tab, angles, irr_meas, fixed, sensor_in, apriori)

    if ~isempty(apriori)
        tab(strcmp(tab.variable,'SMC'), 'value') = apriori.SMCm;
        tab(strcmp(tab.variable,'LAI'), 'value') = apriori.LAIm;
        tab(strcmp(tab.variable,'FVC'), 'value') = apriori.FVCm;

        tab(strcmp(tab.variable,'LAI'), 'uncertainty') = apriori.LAIs;
        tab(strcmp(tab.variable,'FVC'), 'uncertainty') = apriori.FVCs;
    end
    
    if sensor_in.update_Kb == 0
        %% measurement error (see Rodgers 2000, Sect 3)
        epsilon = measurement.std; % instrumental noise
        epsilon(isnan(epsilon)) = Inf;
        iparams = tab.include;

        tab.include = ~tab.include;
        tab.include(contains(tab.variable, ["SIF","V2Z"])) = 0;
        Sb = tab.uncertainty(tab.include);

        % here modification of parameters should occur
        tab = helpers.modify_tab_parameters(tab);
        
        % save prior (x0) values
        tab.x0 = tab.value;
    
        Kb = helpers.clc_jacobian(-999, measurement, tab, angles, irr_meas, fixed, sensor_in);
        measurement.std = diag((diag(epsilon.^2) + Kb*diag(Sb.^2)*Kb.').^0.5); % noise + RTMo parameter error
        tab.include = iparams;
    else
        %% here modification of parameters should occur
        tab = helpers.modify_tab_parameters(tab);
        
        %% save prior (x0) values
        tab.x0 = tab.value;

        iparams = tab.include;
    end
    
    %% initial parameter values, and boundary boxes
    params0 = tab.value(iparams);
    lb = tab.lower(iparams);
    ub = tab.upper(iparams);
    
    stoptol = 1E-6;  % we recommend e-6
    opt = optimset('MaxIter', 30, 'TolFun', stoptol, ...
                   'DiffMinChange', 1E-4); % works for float32 input
                    % 'Display', 'iter');
                   
    
    %% function minimization
    f = @(params)COST_4SAIL_common(params, measurement,  tab, angles, irr_meas, fixed, sensor_in);

    if any(tab.include)  % analogy of any(include == 1)
        tic
        [paramsout,~,~,exitflag,output,~,Jac]= lsqnonlin(f, params0, lb, ub, opt);
        toc
    else % skip minimization and get resuls of RTMo_lite run with initial  parameters (param0)
        paramsout = params0;
    end
    
    %% best-fitting parameters
    tab.value(tab.include) = paramsout;
    results.parameters = helpers.demodify_parameters(tab.value, tab.variable);

    %% best-fittiing spectra
    [er, rad, reflSAIL, rmse, soil, fluo, meas_std] = f(paramsout);
    
    results.rmse = rmse;
    results.refl_mod = reflSAIL;
    results.sif = fluo.SIF;
    results.sif_norm = fluo.SIFnorm;
    results.soil_mod = soil.refl_in_meas;
    results.exitflag = exitflag;
    results.std_refl = meas_std;
    
end
