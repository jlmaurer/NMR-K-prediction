function [b, sig, likes, accept_rat] = mcmc_nmr_bsig (Dk, T2ML, z, n, Niter, stepsize, figureson)

if nargin < 5
    if nargin < 4
        Niter = 1e6; 
    end
    stepsize = 0.8; 
end

    kk = log10(Dk); 
    
    % set bounds on each parameter
    maxlogb = inf; 
    minlogb = -inf; 
    minsig = 1e-5; 
    maxsig = inf; 

    bounds = [minlogb, maxlogb; minsig, maxsig]; 
    k = 100; 
    bcut = 100; 
    x_init = [-2, 1]'; 

    % determine which equation to use: 
    % ip = 0        % no T2B or porosity, log space

    % run mcmc
    [xhats, likes, lkpreds, accept_rat] = mcmc(Niter, stepsize, @NMRfun2, x_init, k, ...
        bounds, [], [], bcut, [],Dk, T2ML, n);

    b = xhats(1,:); 
    sig = xhats(2,:); 
    
    %% Stats
    [mm, nn] = size(xhats); 
    for i = 1:mm
        [h, x] = hist(xhats(i, :), sqrt(nn));
        ind = h == max(h); 
        if sum(ind) ~= 0
            ind = find(ind ~=0, 1, 'first');
        end
        bestps(i) = x(ind); 
    end
    meanvals = mean(xhats, 2);
    meanvals(3) = 10^meanvals(1);
    sdvals = std(xhats, [], 2); 
    p = prctile(xhats, [2.5, 97.5], 2); 
    pl= p(:,1)'; 
    pu = p(:,2)'; 
    pl(1) = 10^pl(1); 
    pu(1) = 10^pu(1); 

    varnames = {'log10_b', 'sig_d'};
    obsnames = {'MAP_values', 'Mean_values', 'SE', 'CI_Upper', 'CI_Lower'}; 
    parammat = [[10^bestps(1), bestps(2)]; [meanvals(2:3)]'; [sdvals(:)]'; pu; pl]; 
    SumP = dataset(parammat(:,1), parammat(:,2), 'VarNames', varnames, 'ObsNames', obsnames);

    %% Plotting
    if figureson == 1
        names = {'log_{10}(b)', '\sigma'}; 
        graph_correlations(xhats,1, names, 1, 0)
        
        plotKwithz2(kk, z, lkpreds)
    end
end