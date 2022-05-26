function fixed = read_fixed_input()

    %% define spectral regions
    spectral.wlP = 400:1:2400;                      % PROSPECT data range
    spectral.wlE = 400:1:750;                       % excitation in E-F matrix
    spectral.wlF = 640:1:850;                       % chlorophyll fluorescence in E-F matrix

    %% fixed input for SIF simulation with PCA
    PCflu = xlsread(fullfile('..', 'input', 'PC_flu.xlsx'));
    pcf   = PCflu(2:end, 2:5);

    %% fixed input for FLUSPECT and BSM
    optipar = load(fullfile('..', 'input', 'fluspect_data', 'Optipar2017_ProspectD'));
%     optipar.optipar.Ks = table2array(readtable(fullfile('..', 'input', 'fluspect_data', 'Ks_Pacheco-Labrador_etal_2021.txt')));
    
    %% collect
    fixed.spectral = spectral;
    fixed.pcf = pcf;
    fixed.optipar = optipar.optipar;
    fixed.srf_sensors = {'OLCI', 'MSI', 'altum_multispec', 'SLSTR', 'Synergy', 'Synergy_official', 'MODIS_Aqua'};
end
