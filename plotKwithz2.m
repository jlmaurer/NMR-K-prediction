function plotKwithz2 (kk, z, lkpred)

    ydat2 = lkpred(:);
    xdat2 = repmat(z, length(lkpred), 1); 
    m = length(ydat2); 
    
    % reduce plotted data
    outind = 1;
    ind = 1:outind:m; 

    figure; hold on
    plot(ydat2(ind), xdat2(ind), 'color', [0.7 0.7 0.7])
    dscatter(ydat2(ind), xdat2(ind), 'msize', 50)
    plot(kk, z, '-*k')
    [m, n] = size(lkpred); 
    if m > n
        lkpred = lkpred'; 
    end
    p = prctile(lkpred,[2.5, 97.5],2);
    p1 = p(:,1);
    p2 = p(:,2); 
    plot(p1,z, '-.k', p2,z, '-.k')
    xlabel('log_{10}(K)', 'FontSize', 12)
    ylabel('Depth (m)', 'FontSize', 12)
    str = ['Predicted and Measured Permeability vs. Depth'];
    title(str, 'FontSize', 14)
    axis ij
end