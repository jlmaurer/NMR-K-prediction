%% Script to estimate m and sigma at each site

% use each individual well
all_names = {'A1', 'C1', 'dpnmr_leque_east', 'dpnmr_leque_west', ...
    'dpnmr_larned_east', 'dpnmr_larned_west', 'dpnmr_larned_lwph'}; 

% use each site
% all_names = {'dpnmr_larned_all', 'dpnmr_leque_all', 'gems_all'}; 

% use all data together
% all_names = {'all_data'}; 
figureson = 0; 

% [b_mcmc, n_mcmc, m_mcmc, sig_mcmc] = deal(zeros(
for k = 1:length(all_names)
    name = all_names{k}; 
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(name); 

    %% MCMC for solution to various parameters
    Niter= 1e6; 
    stepsize = 0.8; 
    n = 2; 

    % MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
    [paramhats, ~, ~, ~, accept_rat] = mcmc_nmr_full(Dk, T2ML, phi, ...
        z, Niter, stepsize, figureson);
    b_mcmc(:,k) = paramhats(2,:); 
    n_mcmc(:,k) = paramhats(3,:);
    m_mcmc(:,k) = paramhats(4,:);
    sig_mcmc(:,k) = paramhats(5,:);
end

%% Plotting
dnames = {'A1', 'C1', 'Leque East', 'Leque West','Larned East', ...
    'Larned West', 'Larned C'}; 

figure; 
subplot(221)
for k = 1:length(all_names)
    histogram(n_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
xlabel('n', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)

subplot(222)
for k = 1:length(all_names)
    histogram(m_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
xlabel('m', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
xlim([-2.1, 5])
set(gca, 'FontSize', 12)

subplot(223)
for k = 1:length(all_names)
    histogram(sig_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
xlim([0.1, 1.8])
xlabel('\sigma', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)
legend(dnames, 'FontSize', 14)

