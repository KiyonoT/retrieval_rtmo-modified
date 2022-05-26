function J = clc_jacobian(fx, measurement, tab, angles, irr_meas, fixed, sensor)
    x = tab.value(tab.include);
    n = length(x); 
    if fx == -999
        fx = REFL_4SAIL(x, measurement, tab, angles, irr_meas, fixed, sensor);
    end
    step = 1e-6; 
    J = zeros(length(fx), n);

    for k=1:n
        xstep = x;
        xstep(k)=x(k) + step;
        fxstep = REFL_4SAIL(xstep, measurement, tab, angles, irr_meas, fixed, sensor); 
        J(:,k)= (fxstep - fx) ./ step;
    end
end

function refl = REFL_4SAIL(p, measurement, tab, angles, irr_meas, fixed, sensor)
    % REFL_4SAIL
    % RETURNS
    %   refl    [measurement.wl]    reflectance with user's input in measurements wavelength
    
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
        [refl, soil_refl, ~] = to_sensor.rtmo2srf_brf(rad, SIF, soil.refl, spectral.wlP, irr_meas, sensor);
%         [refl, soil_refl, fluo] = to_sensor.rtmo2srf_refl(rad, SIF, soil.refl, spectral.wlP, irr_meas, sensor);
%         [refl, soil_refl, fluo] = to_sensor.rtmo2srf_radiance(rad, SIF, soil.refl, spectral.wlP, irr_meas, sensor);
    else
        [refl, soil_refl, ~] = to_sensor.rtmo2measured_refl(rad, SIF, soil.refl, spectral.wlP, irr_meas, measurement, sensor.Rin);
    end
    if ~any(strcmp(tab.variable, 'FVC'))
        canopy.FVC = 1;
    end

    refl = refl*canopy.FVC + soil_refl*(1-canopy.FVC);
end
