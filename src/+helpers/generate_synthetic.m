function generate_synthetic(n_spectra, noise_times)
    
    if nargin == 0
        n_spectra = 10;
        noise_times = 0;
    end

    if nargin == 1
        noise_times = 0;
    end

    [measured, tab, angles, irr_meas, fixed, sensor] = helpers.synthetic_input();
    
    outdir = fullfile('..', 'synthetic', [sensor.instrument_name, '_', num2str(n_spectra)]);
%     assert(~exist(outdir, 'dir'), ['directory with synthetic input already exists at %s\n'...
%         'please, rename or delete it and rerun the function'], outdir)
    mkdir(outdir);
    
    %% randomly sample parameters 
    tab.x0 = tab.value;
    iparams = tab.include;
    varnames = tab.variable(iparams);
    lb = tab.lower(iparams);
    ub = tab.upper(iparams);
    
    if sum(iparams) == 1
        params = linspace(lb, ub, n_spectra);
    else
        rng(0, 'twister')  % setting seed == 0
        params = (ub-lb) .* rand(sum(iparams), n_spectra) + lb;
    end
    
    if any(strcmp('LIDFa' , varnames))
        % abs(LIDFa + LIDFb) <= 1
        i_lidfa = strcmp('LIDFa', varnames);
        i_lidfb = strcmp('LIDFb', varnames);
        lidfa = params(i_lidfa, :);
        lidfb = params(i_lidfb, :);
        params(i_lidfa, :) = (lidfa + lidfb) / 2;
        params(i_lidfb, :) = (lidfa - lidfb) / 2;
    end
    
    if verLessThan('matlab', '9.1')  % < 2016b
        varnames_in = '';
    else
        varnames_in = strjoin(varnames, ', ');
    end
    fprintf('Sampled %i parameters: %s\n', length(varnames), varnames_in)
    
    %% run model
    refls = zeros(length(measured.wl), n_spectra);
    tab.value = helpers.modify_parameters(tab.value, tab.variable);  % LAI is complex otherwise
    for i=1:n_spectra
        % disp(i)
        p = params(:, i);
        p = helpers.modify_parameters(p, varnames);
        [er, rad, refl, rmse, soil, fluo] = COST_4SAIL_common(p, measured, tab, angles, ...
                                                              irr_meas, fixed, sensor);
        refls(:, i) = refl;
    end
    
    %% add noise
    % SNR noise (Synergy and OLCI)
%     snr_synergy = csvread("D:\PyCharm_projects\gsa_rtmo_6S\for_paper\retrievability\SNRs\SNRs_S3A.csv", 0, 1);
%     single_noise = refls ./ snr_synergy(1:size(refls, 1));
%     %noise = single_noise * noise_times; constant % of specific noise
%     noise = rand(size(refls)) .* single_noise * noise_times;
    % constant noise up to certain % from measured
    single_noise = (refls * noise_times / 100);
    noise = rand(size(refls)) .* single_noise * 2;  % *2 because mean(rand) = 0.5
    
    %% write data
%     csvwrite(fullfile(outdir, 'synthetic_toa.csv'), refls)
%     csvwrite(fullfile(outdir, 'synthetic_noise.csv'), noise)
    refls = refls + noise;
    
%     validation = [varnames, array2table(params)];
    not_fit_names = tab.variable(~tab.include);
    not_fit_vals = helpers.modify_parameters(tab.value(~tab.include), not_fit_names);
    validation = [[varnames; not_fit_names], array2table([params; repmat(not_fit_vals, 1, size(params, 2))])];

    csvwrite(fullfile(outdir, 'synthetic.csv'), refls)
    csvwrite(fullfile(outdir, 'synthetic_wl.csv'), measured.wl)
    writetable(validation, fullfile(outdir, 'synthetic_val.csv'))
    writetable(tab, fullfile(outdir, 'synthetic_input.csv'))
    
    fid = fopen(fullfile(outdir, 'synthetic_comment.txt'), 'w');
    fprintf(fid, 'sensor %s\n', sensor.instrument_name);
    fprintf(fid, 'angles tts=%.2g, tto=%.2g, psi=%.2g\n', angles.tts, angles.tto, angles.psi);
    fprintf(fid, 'n_spectra=%d\n', n_spectra);
    fclose(fid);
    
    fprintf('saved all in %s\n', outdir)
    
    %% plot results    
    [~, order] = sort(measured.wl);  % for synergy OLCI, SLSTR
    figure()
    plot(measured.wl(order), refls(order, :), 'o-')
    saveas(gcf, fullfile(outdir, [sensor.instrument_name, '_', num2str(n_spectra) '.png']))
%     hold on
%     plot(measured.wl, soil.refl_in_meas, 'x')
%     l = legend([arrayfun(@num2str, params, 'UniformOutput', false), 'soil'], ...
%         'location', 'eastoutside');
%     l.ItemHitFcn = @hitcallback_ex1;

end

function hitcallback_ex1(src,evnt)

if strcmp(evnt.Peer.Visible,'on')
    evnt.Peer.Visible = 'off';
else 
    evnt.Peer.Visible = 'on';
end

end
