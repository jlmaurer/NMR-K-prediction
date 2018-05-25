%% Script for interpolating porosity measurements. 

clear
name = 'all_data'; 
% load data file
[d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(name); 

%  figure; scatter(phi, Dk, [], lt)
% set(gca, 'yscale', 'log')
% colormap(jet)
% colorbar
ind = Dk > 0.0008; 
Dk(ind) = []; 
phi(ind) = []; 
lt(ind) = []; 
T2ML(ind) = []; 
lp(ind) = []; 
kk(ind) =[]; 

% prediction from porosity
phispace = linspace(min(phi), max(phi), 300); 
Kp1 = 1e-2*phispace.^4; 
Kp2 = 3e-4*phispace; 

figure; 
scatter(phi, Dk, [], lt, 'filled')
colormap(jet)
colorbar
hold on
% plot((phispace), Kp1, '-k')
% plot((phispace), Kp2, '--k')
set(gca, 'yscale', 'log')
xlabel('\phi')
ylabel('log_{10}(K)')


% figure; 
% scatter3(phi, lt, kk, [],phi, 'filled')
% xlabel('\phi')
% ylabel('log_{10}(T_{2ML})')
% zlabel('log_{10}(K)')
% % Grid to interpolate
% ng = 500; 
% x = linspace(min(lt), max(lt), ng); 
% y = linspace(min(kk), max(kk), ng); 
% [X, Y] = meshgrid(x,y); 
% 
% PHI = griddata(lt, kk, phi, X(:), Y(:), 'natural'); 
% 
% figure; 
% scatter3(lt, kk, phi, 'filled')
% hold on
% scatter3(X(:), Y(:), PHI(:), [], PHI(:), '.')
