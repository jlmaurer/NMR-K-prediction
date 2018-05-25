function [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full (Dk, T2ML, phi, z, Niter, stepsize, figureson)
kk = log10(Dk); 

% parameters to search over: T2B, b, n, m
% set bounds on each parameter
maxT2B = 15; 
minT2B = 0; 
maxlogb = inf; 
minlogb = -inf; 
maxn = 4; 
minn = 0; 
maxm = 5; 
minm = -2; 
minsig = 1e-5; 
maxsig = inf; 

bounds = [minT2B, maxT2B; minlogb, maxlogb; minn, maxn; ...
    minm, maxm; minsig, maxsig]; 
k = 100; 
bcut = 100; 
x_init = [2.5, -2, 2, 0, 1]'; 

% determine which equation to use: 
% ip = 0        % no T2B or porosity, log space
% ip = 1        % T2B and porosity, log space
% ip = 2        % T2B, no porosity, log space
% ip = 3        % no T2B or porosity, real space - not stabile

ip = 1;

% run mcmc
[paramhats, likes, kpreds, accept_rat] = mcmc(Niter, stepsize, @NMRfun, x_init, k, ...
    bounds, [], [], bcut, [],Dk, phi, T2ML, ip);

if figureson == 1
    names = {'T_{2B}', 'log_{10}(b)', 'n', 'm', '\sigma'}; 
    graph_correlations(paramhats,1, names, 1)

% %%%%%%%%%%%%%%%% Can also plot all models with depth
%     [ksort, ik] = sort(kk); 
%     xdat = repmat(ksort, length(kpreds), 1); 
%     ydat = kpreds(ik, :);
%     figure; hold on
%     % plot(ksort, ydat, 'color', [0.7, 0.7, 0.7])
%     ydat = ydat(:); 
%     dscatter(xdat, ydat, 'msize', 60)
%     plot(kk, kk, '-*k')
%     p = prctile(kpreds,[2.5, 97.5],2);
%     p1 = p(ik,1); 
%     p2 = p(ik,2); 
%     plot(ksort, p1, '-.k', ksort, p2, '-.k')
%     xlabel('Measured log_{10}(k)', 'FontSize', 12)
%     ylabel('Estimated log_{10}(k)', 'FontSize', 12)
%     str = ['Predicted vs. Measured Permeability'];
%     title(str, 'FontSize', 14)
% 
%     xdat2 = repmat(z, length(kpreds), 1); 
%     figure; hold on
%     % plot(kpreds,z, 'color', [0.7 0.7 0.7])
%     ydat2 = kpreds(:); 
%     dscatter(ydat2(1:11:end), xdat2(1:11:end), 'msize', 60)
%     plot(kk,z, '-*k')
%     p1 = p(:,1);
%     p2 = p(:,2); 
%     plot(p1,z, '-.k', p2,z, '-.k')
%     axis ij
%     xlabel('log_{10}(K)', 'FontSize', 12)
%     ylabel('Depth (m)', 'FontSize', 12)
%     str = ['Predicted and Measured Permeability vs. Depth'];
%     title(str, 'FontSize', 14)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% find MAP solution
[mm, nn] = size(paramhats); 
bestps = zeros(mm, 1); 
lci = zeros(mm, 1); 
hci = zeros(mm,1); 
for i = 1:mm
    [h, x] = hist(paramhats(i, :), sqrt(nn));
    ind = h == max(h);
    if sum(ind) ~= 1
        ind = find(ind ~= 0, 1, 'first');
    end
    bestps(i) = x(ind); 
    p = prctile(paramhats(i,:), [2.5, 97.5]); 
    lci(i) = p(1); 
    hci(i) = p(2); 
end

meanvals = mean(paramhats, 2); 
stdvals = std(paramhats,[], 2); 
minvals = min(paramhats,[], 2); 
maxvals = max(paramhats,[], 2); 

Varnames = {'T2B', 'logb', 'b', 'n', 'm', 'Sigma'}; 
obsnames= {'MAP_values', 'Mean_values', 'SE', 'Minimum', 'Maximum', 'CI_Lower', 'CI_Upper'}; 
datamat = [bestps'; meanvals'; stdvals'; minvals'; maxvals'; lci'; hci'];
hps = dataset(datamat(:,1), datamat(:,2), 10.^datamat(:,2), datamat(:,3), datamat(:,4), datamat(:,5), 'VarNames', Varnames, 'ObsNames', obsnames); 

%% Export file if wanted
% filename = ['Summary_Parameters_for_' name '.txt'];
% export(hps, 'file', filename)

end