all_names = {'A1','C1', 'dpnmr_leque_east', 'dpnmr_leque_west', 'dpnmr_larned_east', ...
    'dpnmr_larned_west', 'dpnmr_larned_lwph', 'gems_all', 'dpnmr_leque_all','dpnmr_larned_all',  'all_data'}; 


C = {'-*k', '-*r', '-.b', '-.g', '-oc', '-om', '-or', '-*b', '-.b', '-ob', '--k'}; 
figure; 
hold on

for k = 1:length(all_names)

    clearvars -except all_names k rnorm C
    
    name = all_names{k}; 

    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(name); 
    logSumEch = log10(SumEch); 

    %% Test which SOE works best
    % 1 - SDR
    % 2 - SOE
    % 3 - SOE - 3 seconds (extrapolated)
    % 4 - SOE - time-weighted mean
    % 5 - SOE - 3 second time-weighted mean

    vars = [lt, log10(SumEch), log10(SumEch_3s), log10(SumEch_twm), log10(SumEch_twm_3s)]; 
    numvars = size(vars, 2);
    rnorm = []; 
    for i = 1:numvars
        var = vars(:,i); 
        r = kk - (var\kk)*var; 
        rnorm(i) = norm(r)./length(r); 
    end
    
    plot([1:numvars], rnorm, C{k})

end
set(gca, 'xtick', [1:5])
labels = {'T_{2ML}', 'SOE', 'SOE - 3s', 'SOE - TWM', 'SOE - 3s, TWM'}; 
set(gca, 'xticklabel', labels)
xlim([0.5, 5.5])
ax = gca;
ax.XTickLabelRotation = 45;

names = {'GEMS2 - A1', 'GEMS2 - C1', 'Leque East', 'Leque West', 'Larned East', 'Larned West', 'Larned C', 'GEMS2 - all', 'Leque - all', 'Larned - all', 'All Data'}; 
legend(names')


%% Now look at SOE
Nboot =  200; % number of bootstrap samples

for k = 1:length(all_names)

    name = all_names{k}; 

    % load data file
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(name); 
    logSumEch = log10(SumEch); 
    logSumEch_3s = log10(SumEch_3s); 
    logSumEch_twm = log10(SumEch_twm); 
    logSumEch_twm_3s = log10(SumEch_twm_3s); 
    %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
    lt = logSumEch_3s; 
    T2ML = SumEch; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [b_boot, n_boot] = bootstrap_fun_mb([lt, kk], Nboot);    % n can vary
    
    meanb(k) = mean(b_boot); 
    meann(k) = mean(n_boot); 

end

figure; 
plot([1:k], meann)
