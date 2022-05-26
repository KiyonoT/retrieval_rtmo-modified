function [er, rad, refl, rmse, soil, fluo, meas_std] = COST_4SAIL_common(p, measurement, tab, angles, ...
                                                               irr_meas, fixed, sensor)
    % COST_4SAIL_my
    % RETURNS
    %   er                          difference between modeled and measured reflectance + prior weight
    %   rad     [wlP]               reflectance factors of RTMo in PROSPECT wavelength
    %   refl    [measurement.wl]    reflectance with user's input in measurements wavelength
    %   rmse    [wlmin:wlmax]       RMSE of modeled-measured in fit wavelength
    %   soil    [wlP, measurement.wl] soil reflectance in PROSPECT wavelength and mwasurements wl
    
    %% fixed parameters
    
    spectral = fixed.spectral;
    optipar = fixed.optipar;
    pcf = fixed.pcf;
    
    %% create input structs from table
    tab.value(tab.include) = p;
    tab.value = helpers.demodify_parameters(tab.value, tab.variable);

    soilpar = io.table_to_struct(tab, 'soil');
    canopy = io.table_to_struct(tab, 'canopy');
    leafbio = io.table_to_struct(tab, 'leafbio');
    wpcf = io.table_to_struct(tab, 'sif');

    %% leaf reflectance - Fluspect
    leafbio.fqe(2)  = 0.02;                     % quantum yield
    leafbio.fqe(1)  = 0.02 / 5;
    if ~any(strcmp(tab.variable, 'V2Z'))
        leafbio.V2Z = 0;
    end

%     leafopt = models.fluspect_B_CX_PSI_PSII_combined(spectral, leafbio, optipar);
    leafopt = models.fluspect_lite(spectral, leafbio, optipar);

    %% soil reflectance - BSM
    soilemp.SMC   = 25;        % empirical parameter (fixed) [soil moisture content]
    soilemp.film  = 0.015;     % empirical parameter (fixed) [water film optical thickness]
    % soilspec.wl  = optipar.wl;  % in optipar range
    soilspec.GSV  = optipar.GSV;
    soilspec.kw   = optipar.Kw;
    soilspec.nw   = optipar.nw;
    
    if isfield(measurement, 'soil_refl')
        soil.refl = measurement.soil_refl;
%         soil.refl = interp1(measurement.wl, measurement.soil_refl,
%         spectral.wlP, 'splines', 1E-4); % 'spline' since M2020
    else
        soil.refl = models.BSM(soilpar, soilspec, soilemp);
    end

    %% canopy reflectance factors - RTMo
    canopy.nlayers  = 60;
    nl              = canopy.nlayers;
    canopy.x        = (-1/nl : -1/nl : -1)';         % a column vector
    canopy.xl       = [0; canopy.x];                 % add top level
    canopy.nlincl   = 13;
    canopy.nlazi    = 36;
    canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
    canopy.lazitab  = ( 5:10:355 );           % a row
    canopy.hot      = sensor.hot;
    canopy.lidf     = equations.leafangles(canopy.LIDFa, canopy.LIDFb); 

    rad   = models.RTMo_lite(soil, leafopt, canopy, angles);
%     rad   = models.RTMo_lite_c(soil, leafopt, canopy, angles);
    
    %% canopy fluorescence from PCA, in W m-2 sr-1
    SIF = zeros(length(spectral.wlP),1);
    SIF(640-399:850-399)   = pcf * cell2mat(struct2cell(wpcf));
    rad.SIF = SIF(640-399:850-399);

    %% canopy reflectance in measurements wl
    if isfield(sensor, 'srf')
        [refl, soil_refl, fluo] = to_sensor.rtmo2srf_brf(rad, SIF, soil.refl, spectral.wlP, irr_meas, sensor);
%         [refl, soil_refl, fluo] = to_sensor.rtmo2srf_refl(rad, SIF, soil.refl, spectral.wlP, irr_meas, sensor);
%         [refl, soil_refl, fluo] = to_sensor.rtmo2srf_radiance(rad, SIF, soil.refl, spectral.wlP, irr_meas, sensor);
        if isfield(measurement, 'xa')
            y = refl ./ (1 - refl .* measurement.xc);
            toa = (y + measurement.xb) ./ measurement.xa;
            er1 = toa - measurement.rad;
        elseif any(strcmp(tab.variable, 'FVC'))
            er1 = (refl*canopy.FVC + soil_refl*(1-canopy.FVC)) - measurement.refl;
        else
            er1 = refl - measurement.refl;
        end
    else
        [refl, soil_refl, fluo] = to_sensor.rtmo2measured_refl(rad, SIF, soil.refl, spectral.wlP, irr_meas, measurement, sensor.Rin);
        I = measurement.i_fit;
        er1 = refl(I) - measurement.refl(I);
    end
    
    soil.refl_in_meas  = soil_refl;
    
    %% extra weight from prior information
    prior.Apm = tab.x0(tab.include);
    prior.Apm(strcmp(tab.variable(tab.include),'Cca')) = 0.25*leafbio.Cab; % Note: 1 + 0.22*Cab (LOPEX + ANGERS); 2 + 0.19*Cab (Sims & Gamon 2002)
    prior.Aps = tab.uncertainty(tab.include);
    er2 = (p - prior.Apm) ./ prior.Aps; 
    
    %% calculate the difference between measured and modeled data
    if sensor.update_Kb == 0
        meas_std = measurement.std;
    else
        epsilon = measurement.std; % instrumental noise
        epsilon(isnan(epsilon)) = Inf;
        
        % add model error (see Rodgers 2000, Sect 3)
        tab.include = ~tab.include;
        tab.include(contains(tab.variable, ["SIF","V2Z"])) = 0;
        Sb = tab.uncertainty(tab.include);
        tab = helpers.modify_tab_parameters(tab);
        Kb = helpers.clc_jacobian(refl*canopy.FVC+soil_refl*(1-canopy.FVC), ...
                                  measurement, tab, angles, irr_meas, fixed, sensor);
        meas_std = diag((diag(epsilon.^2) + Kb*diag(Sb.^2)*Kb.').^0.5); % noise + RTMo parameter error
    end
    
    er = er1(~isnan(er1)) ./ meas_std(~isnan(er1));
    er1 = er1(~isnan(er1));
    
    %% total error
    er = [er ; er2];

    rmse = sqrt((er1' * er1) ./ numel(er1));

end
