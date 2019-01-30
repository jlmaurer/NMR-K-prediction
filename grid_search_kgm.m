function [tauspace, rhospace, r] = grid_search_kgm2(lT2, lphi, kk)
% this function does a grid search over specified parameter values in the
% SDR equation

% put in linear space
T2 = 10.^lT2; 
phi = lphi; 
K = kk; 

% temp-dependent parameters
Tb= @(T) 3.3 + .044.*(T - 35);                  % Bulk relaxation time constant (s)
eta = @(T) 2.414e-5*(10.^((247.8)./(T + 133.15)));  % viscosity (Pa s)
aa = 1.0413; bb = 0.039828; cc = 0.00040318;            % diffusivity constants
D = @(T) aa + bb*T+ cc*T.^2;   % diffusivity (m^2/s)
rho_h2o = @(T) 1000*(1 - (((T + 288.94)/(508929 ...
    *(T+68.12)))*((T - 3.98).^2)));             % density of water (kg/m^3)

% set bounds on A = g/(8*tau^2*rho^2), B = rho^2
max_tau = 0; 
min_tau = -2; 

max_rho = 0; 
min_rho = -6.6; 

nm = 201; 
tauspace = linspace(min_tau, max_tau, nm); 
rhospace = logspace(min_rho, max_rho, nm); 

% KGM equation
% lK_kgm = @(A,B, T, T2, lphi) (A) + lphi ...
%     + 2*log10( - (D(T)/B) + sqrt((D(T)/B)^2 + ((4*D(T)*Tb(T)*T2)./(Tb(T) - T2)))); 
lK_kgm = @(A, B, T, T2, phi) KGMfun([A, B], [T2]); 

T = 20; 
% run grid search
[r] = deal(zeros(length(tauspace), length(rhospace)));  
for mloop = 1:length(rhospace)
    btest = rhospace(mloop);    
    for bloop = 1:length(tauspace)    
        atest = tauspace(bloop);  
        lKp = lK_kgm(atest, btest, T, T2, phi);
        r(mloop, bloop) = norm(lKp - kk); 
    end
end

% % compute best-fitting values
% ind1 = Bspace == 1; 
% ind2 = Bspace == 2; 
% ind4 = Bspace == 4; 
% ind0 = Bspace == 0; 
% 
% r1 = r(ind1, :) == min(r(ind1,:)); 
% r2 = r(ind2, :) == min(r(ind2,:)); 
% r4 = r(ind4,:) == min(r(ind4,:)); 
% r0 = r(ind0,:) == min(r(ind0,:)); 
% 
% bestb(1)= log10(Aspace(r1)); 
% bestb(2)= log10(Aspace(r2)); 
% bestb(3)= log10(Aspace(r4)); 
% bestb(4)= log10(Aspace(r0)); 

s = size(r); 
bestp = min(r(:)); 

indp = find(r(:) == bestp); 
[ii, jj] = ind2sub(s, indp);
bestA = tauspace(ii); 
bestB = rhospace(jj); 
lkpred = lK_kgm(bestA, bestB, T, T2, lphi); 

figure; scatter(kk, lkpred), axis equal, hold on, refline(1,0)
xlabel('Observerd K')
ylabel('Predicted K')

taus = 1./sqrt(10.^(tauspace)); 
% Plot grid results
figure; 
% subplot(4, 4, [1:3,5:7,9:11, 13:15])
hold on
imagesc(taus, log10(rhospace), log10(r/length(phi)))
plot(taus(jj), log10(rhospace(ii)), '*w')
caxis([-1.2,0])
% plot( bestb, [1, 2, 4, 0], 'sw')
colorbar('Ticks', [-1.2:.2:0])
% colorbar
colormap(jet)
xlim([min((taus)), max((taus))])
ylim([min(log10(rhospace)), max(log10(rhospace))])
ylabel('log_{10} (\rho)')
xlabel('\tau')
% str = ['(a)  Misfit - best value pair: m = ', num2str(Bspace(ii)) ', log(b) = ',num2str(log10(Aspace(jj)))];
% title(str)

% subplot(4,4,16)
% hold on
% plot(log10(Aspace), r(ind0,:))
% title('(e)  m = 0')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% subplot(4,4,12)
% hold on
% plot(log10(Aspace), r(ind1,:))
% title('(d)  m = 1')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% subplot(4,4,8)
% hold on
% plot(log10(Aspace), r(ind2,:))
% title('(c)  m = 2')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% subplot(4,4,4)
% hold on
% plot(log10(Aspace), r(ind4,:))
% title('(b)  m = 4')
% xlabel('log_{10}(b)')
% ylabel('Residual')
% 
% set(gcf,'OuterPosition',[100 100 1000 800], 'PaperPositionMode','auto','color',[1 1 1])

% % output best b values
% bs = [Aspace(ii); bestb(:)]; 
% ms = [Bspace(jj); [0,1, 2, 4]']; 
end