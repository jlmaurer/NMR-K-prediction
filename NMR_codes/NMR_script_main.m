%% Main script for analysis of NMR data
% This script analyzes a set of NMR data using the grid search, bootstrap,
% stepwise least squares regression, and MCMC. 
%
% text files that are called here were created using the script
% 'NMR_textfile_build_script' callin 'build_nmr_txt_file' function. This
% function gets the data found in bstrp_dat and adds in the extra
% Sum-of-Echoes data found in the folders inside bstrp_dat_VC, and outputs
% the new text file to bstrp_dat_VC_processed. These test files are what
% are called here. NOTE: can also revert back to calling the files directly
% from bstrp_dat by editing the load data command to 'loadnmrdata(name)'. 
% Note that I leave the permeability data in the original Vista Clara form,
% to be converted to correct units by loadnmrdata2.m. 

%% Data Details
% data should be a text file with columns given by depth,
% T2ML, porosity, permeability, and Vista Clara units
%
%% Options for Data files
%   'A1'                    % this is a gems2 location
%   'C1'                    % this is a gems2 location
%   'gems_all'              % both A1 and C1 - using Rosemary's names
%   'dpnmr_larned_east'
%   'dpnmr_larned_west'
%   'dpnmr_larned_lwph'
%   'dpnmr_larned_all'
%   'dpnmr_leque_east'
%   'dpnmr_leque_west'
%   'dpnmr_leque_all'
%   'all_data'

%% Variable Names
% Dk - Permeability
% T2ML - T2ML
% phi - porosity
% z - depth
% SumEch - raw sum of echoes
% kk - log10(Dk)
% lt - log10(T2ML)
% lp - log10(porosity)

%% Load Data
clear

name = 'A1'; 

% load data file
[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(name); 
logSumEch = log10(SumEch); 
logSumEch_3s = log10(SumEch_3s); 
logSumEch_twm = log10(SumEch_twm); 
logSumEch_twm_3s = log10(SumEch_twm_3s); 
    
% %%%%%%%%% Change T2 variable to Sum of Echoes for the inversions. 
% logT2ML = logSumEch_3s; 
% T2ML = SumEch; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Turn figures on or off
figureson = 1; 
if figureson == 1
    figure; scatter(logT2ML, logK, 'o')
    nums = 1:length(logT2ML);
    text(logT2ML, logK, num2str(nums'))
    xlabel('Log10 T2ML')
    ylabel('Log10 K')
end

%% Test which SOE works best
% 1 - SDR
% 2 - SOE
% 3 - SOE - 3 seconds (extrapolated)
% 4 - SOE - time-weighted mean
% 5 - SOE - 3 second time-weighted mean

vars = [logT2ML, log10(SumEch), log10(SumEch_3s), log10(SumEch_twm), log10(SumEch_twm_3s)]; 
numvars = size(vars, 2);
for i = 1:numvars
    var = vars(:,i); 
    r = logK - (var\logK)*var; 
    rnorm(i) = norm(r); 
end
figure; plot([1:numvars], rnorm)
set(gca, 'xtick', [1:5])
xlim([0.5, 5.5])

%% Basic solving for b for fixed n, m
% Given m and n, we can solve directly for b -> b = log(k) - m*log(phi) -
% n*log(T2ML). This is the 'direct' method.

C = @(m, n, lt, lp) m*lp + n*lt; 
n = 1;
m= 0; 
bdir_n1 = logK - C(m, n, logT2ML, logPhi); 

n = 2; 
bdir_n2 = logK - C(m, n, logT2ML, logPhi); 

if figureson == 1
    figure; 
    subplot(211)
    hist(bdir_n1, floor(sqrt(length(bdir_n1))))
    xlabel('b')
    ylabel('Frequency')
    title('b values estimated directly for fixed n, n = 1')
    subplot(212)
    hist(bdir_n2, floor(sqrt(length(bdir_n2))))
    xlabel('b')
    ylabel('Frequency')
    title('b values estimated directly for fixed n, n = 2')
end

%% Linear fit inversion, solving for m, n, and b
% Here I perform a stepwise least sqaures regression (machine learning)
% algorithm, as well as two ordinary least squares regressions. The
% stepwise algorithm takes only those variables it finds improves the fit
% to the model. 

LM1 = fitlm([logPhi, logT2ML], logK, 'linear');
LM2 = fitlm(logT2ML, logK, 'linear'); 
LM3 = stepwiselm([logPhi, logT2ML], logK, 'interactions'); 

if figureson ==1
    figure; 
    subplot(211)
    plotResiduals(LM1)
    title('Histogram of Residuals using T_{2ML} and \phi')
    subplot(212)
    plotResiduals(LM2)
    title('Histogram of Residuals using T_{2ML} only')
end

bs = logspace(-6, 0.3, 1000); 
Y = lognpdf(bs, LM2.Coefficients.Estimate(1), LM2.Coefficients.SE(1)); 
figure;
plot(bs, Y, 'm', 'Linewidth', 2)
xlabel('b')
ylabel('Predicted K')

%% Grid search to test for dependence on porosity
% systematically vary m and b to find the best-fitting paramater region.
% Produces a color plot showing parameter value-pair residuals. 

if figureson == 1
    n = 2;      % fix n for the grid search, allow m to vary
    [bs, ms, r] = grid_search(logT2ML, logPhi, logK, n);
end

%% Bootstrap
% Randomly sample data with replacement, solve the subset for the
% best-fitting parameter values, and repeat many times. 

% m assumed 0; 
n = 2; 
m = 1; 
Nboot =  2000; % number of bootstrap samples

% Takes [log10(T2ML), log10(K)] or [log10(T2ML), log10(phi), log10(K)] as a
% single matrix
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt,lp, kk], Nboot);         % m, n can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt, lp, kk], Nboot, n);        % m can vary
[b_boot, n_boot, m_boot] = bootstrap_fun_mb([logT2ML, logK], Nboot);    % n can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt,lp, kk], Nboot, n, m);   % m, n fixed

if figureson ==1 
    bs = log10(b_boot); 
    graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
end
meanb = mean(b_boot)
sortb = sort(b_boot); 
blo = sortb(50)
bhi = sortb(1950)

%% MCMC for solution to various parameters
% Markov Chain Monte Carlo using Metropolis-Hastings algorithm. Assumes
% Bayes' theorem. Parameters Niter (number of iterations) and stepsize may
% need to be adjusted for good covergence, although with this data set the
% values shown perform reasonably well.
Niter= 1e6; 
stepsize = 0.8; 

%%%%%%%%%%%%%%%%%%%% MCMC hps gives summary parameters (average, CI, etc.) Pick
% either this one or simpler MCMC algorithm to compare with other results.
% NOTE: if the first MCMC algorthm is used, remember when comparing results
% between different methods that the others use fixed values for m and n,
% while this allows both to vary. This will inevitably affect the range in
% b obtained. 

%%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
[paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
    z, Niter, stepsize, figureson);
T2B_mcmc  = paramhats(1,:);
b_mcmc = paramhats(2,:); 
n_mcmc = paramhats(3,:);
m_mcmc = paramhats(4,:);
sig_mcmc = paramhats(5,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%% MCMC for b and data error only. n is fixed and m is assumed zero. 
n = 2;  % fixed n
[b_mcmc, sig_mcmc, likes, accept_rat] = mcmc_nmr_bsig (K, T2ML, z, n, ...
    Niter, stepsize, figureson); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if figureson == 1
    graph_correlations(paramhats, 2, {'T_B', 'log_{10}(b)', 'n', 'm', '\sigma'}, 0, 0)
    all_lKpreds = zeros(length(T2ML), length(b_mcmc)); 
    for k = 1:length(b_mcmc)
        bbm = [b_mcmc(k), sig_mcmc(k)]; 
        [~,all_lKpreds(:,k)] = NMRfun2(bbm, K, T2ML, n); 
    end
    figure; 
    hold on
    for k = 1:size(all_lKpreds,1)
        dpk = all_lKpreds(k,:); 
        [h,x] = hist(dpk, 40); 
        t2m = repmat(log10(T2ML(k)), 1, length(x)); 
        scatter(t2m, x, [], h); 
    end
    scatter(log10(T2ML), logK, 'ok', 'filled')
end

%% Compare all results
% create a matrix to plot and compare all results
mublm = LM2.Coefficients.Estimate(1); 
seblm = LM2.Coefficients.SE(1); 
bslm = mublm+ seblm*randn([length(b_boot), 1]); 
X = zeros(length(b_mcmc), 4); 
X(:,1) = 10.^b_mcmc(:); 
X(1:length(b_boot),2) = b_boot(:); 
X(1:length(bdir_n2),3) = 10.^bdir_n2; 
X(1:length(bslm),4) = 10.^bslm; 
zind = X == 0; 
X(zind) = nan; 

%% Save results
savealldata = 0;    % if =1 save results
if savealldata ==1 
    str = ['Final_full_analysis_for_' name]; 
    save(str)
end

%% Plot final results
% need to trim distributions
names = {'MCMC', 'Bootstrap', 'Direct', 'Analytical SE'};
figure; hold on
boxplot(X, names, 'colors', 'rgbk')
ylabel('b')
