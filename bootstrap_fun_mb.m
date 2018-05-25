function [bs, ns, ms, logb] = bootstrap_fun_mb (data, Nboot)

m = 0; 
% call function, depending on whether or not m and n are specified
bhat = bootstrp(Nboot, @(x) bootstrap_nmr_called_m(x, m), data);    % n can vary
    
[~, bbb] = size(bhat); 
if bbb > 3
    bhat = bhat'; 
    [~, bbb] = size(bhat); 
end
if bbb == 3
    bs = bhat(:,1); 
    ns = bhat(:,2);
    ms = bhat(:,3); 
elseif bbb == 2
    bs = bhat(:,1); 
    ns = bhat(:,2);
    ms = []; 
else
    bs = bhat;
    ns = []; 
    ms = []; 
end

% transform back into log space
logb = bs; 
bs = 10.^logb;
end