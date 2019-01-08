function graph_correlations( corr_mat, dim , names, contourson, fiton)
%GRAPH_CORRELATIONS This functions plots variables and shows both the
%marginal distributions along the diagonal and the correlations between
%each variable on the off-diagonals. SHOULD NOT USE FOR MORE THAN ABOUT 10
%VARIABLES.

if nargin ==2
    names = dim; 
elseif nargin == 4
    fiton = 0; 
end
[a, b] = size(corr_mat); 
if a < b
    corr_mat = corr_mat'; 
end
[mrows, ind] = size(corr_mat); 
if ind > 15
    error('Too many variables, chart will not look good')
end
nbins = ceil(sqrt(mrows)); 
bmat = zeros(ind); 

figure; 
for i = 1:ind
    for j = 1:ind
        if i <= j
            h=subplot(ind, ind, sub2ind(size(bmat), j, i));
            if i == j
                histogram(corr_mat(:,i), nbins, 'Normalization', 'pdf')
                if fiton == 1
                    xs = linspace(min(corr_mat(:,i)), max(corr_mat(:,i)), 200); 
                    ys = normpdf(xs, mean(corr_mat(:,i)), std(corr_mat(:,i))); 
                    hold on
                    plot(xs, ys, 'r', 'LineWidth', 2)
                    xlabel(['Fit: \mu = ', num2str(mean(corr_mat(:,i)))])
                end
            else
                dscatter(corr_mat(:,j),corr_mat(:,i), 'marker', 's', 'msize', 4)
                if contourson ==1 
                    hold on
                    dscatter(corr_mat(:,j),corr_mat(:,i), 'plottype', 'contour')
                    hold off
                end
                colormap(jet)
            end
            axis tight
            if i == 1
                title(names{j})
            end
        end
    end
end

end

