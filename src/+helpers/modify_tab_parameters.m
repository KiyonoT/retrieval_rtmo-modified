function tab = modify_tab_parameters(tab)
    %% modify LAI
    transform_lai = @(x) 1 - exp(-0.2 * x);
    i_lai = strcmp(tab.variable, 'LAI');
    if any(strcmp(tab.variable, 'FVC'))
        i_fvc = strcmp(tab.variable, 'FVC');
        tab.value(i_lai) = tab.value(i_lai)/tab.value(i_fvc);
    end

    tab.value(i_lai) = transform_lai(tab.value(i_lai));
    tab.lower(i_lai) = transform_lai(tab.lower(i_lai));
    tab.upper(i_lai) = transform_lai(tab.upper(i_lai));
    tab.uncertainty(i_lai) = transform_lai(tab.uncertainty(i_lai));
    
    %% modify LIDFs
    % abs(LIDFa + LIDFb) <= 1
    % to make it possible LIDFa = LIDFa + LIDFb, LIDFb = LIDFa - LIDFb
    i_lidfa = strcmp(tab.variable, 'LIDFa');
    i_lidfb = strcmp(tab.variable, 'LIDFb');
    lidfa = tab.value(i_lidfa);
    lidfb = tab.value(i_lidfb);
    tab.value(i_lidfa) = lidfa + lidfb;
    tab.value(i_lidfb) = lidfa - lidfb;
    
end
