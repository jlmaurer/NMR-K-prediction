function [bs, ms, r] = grid_search(lT2, lphi, kk, n)
% this function does a grid search over specified parameter values in the
% SDR equation

if nargin < 4
    n = 2; 
end

% set bounds on b, m
max_b = 2; 
min_b = -4; 

max_m = 5; 
min_m = -2; 

nm = 100*(max_m-min_m) +1;
bspace = logspace(min_b, max_b, nm); 
mspace = linspace(min_m, max_m, nm); 

% SDR equation
SDR_f = @(b, m, n, lT2, lphi) log10(b) + m*lphi + n*lT2; 

% run grid search
[r] = deal(zeros(length(bspace), length(mspace)));  
for mloop = 1:length(mspace)
    m = mspace(mloop);    
    for bloop = 1:length(bspace)    
        b = bspace(bloop);  
        Kp = SDR_f(b, m, n, lT2, lphi);
        r(mloop, bloop) = norm(Kp - kk); 
    end
end

% compute best-fitting values
ind1 = mspace == 1; 
ind2 = mspace == 2; 
ind4 = mspace == 4; 
ind0 = mspace == 0; 

r1 = r(ind1, :) == min(r(ind1,:)); 
r2 = r(ind2, :) == min(r(ind2,:)); 
r4 = r(ind4,:) == min(r(ind4,:)); 
r0 = r(ind0,:) == min(r(ind0,:)); 

bestb(1)= log10(bspace(r1)); 
bestb(2)= log10(bspace(r2)); 
bestb(3)= log10(bspace(r4)); 
bestb(4)= log10(bspace(r0)); 

s = size(r); 
bestp = min(r(:)); 
indp = find(r(:) == bestp); 
[ii, jj] = ind2sub(s, indp); 

% Plot grid results
figure; 
subplot(4, 4, [1:3,5:7,9:11, 13:15])
hold on
imagesc(log10(bspace), mspace, r./length(kk))
plot(log10(bspace(jj)), mspace(ii), '*w')
plot( bestb, [1, 2, 4, 0], 'sw')
colorbar
colormap(jet)
xlim([min(log10(bspace)), max(log10(bspace))])
ylim([min(mspace), max(mspace)])
xlabel('log_{10} b')
ylabel('m')
str = ['(a)  Misfit - best value pair: m = ', num2str(mspace(ii)) ', log(b) = ',num2str(log10(bspace(jj)))];
title(str)

subplot(4,4,16)
hold on
plot(log10(bspace), r(ind0,:))
title('(e)  m = 0')
xlabel('log_{10}(b)')
ylabel('Residual')

subplot(4,4,12)
hold on
plot(log10(bspace), r(ind1,:))
title('(d)  m = 1')
xlabel('log_{10}(b)')
ylabel('Residual')

subplot(4,4,8)
hold on
plot(log10(bspace), r(ind2,:))
title('(c)  m = 2')
xlabel('log_{10}(b)')
ylabel('Residual')

subplot(4,4,4)
hold on
plot(log10(bspace), r(ind4,:))
title('(b)  m = 4')
xlabel('log_{10}(b)')
ylabel('Residual')

set(gcf,'OuterPosition',[100 100 1000 800], 'PaperPositionMode','auto','color',[1 1 1])

% output best b values
bs = [bspace(ii); bestb(:)]; 
ms = [mspace(jj); [0,1, 2, 4]']; 
end