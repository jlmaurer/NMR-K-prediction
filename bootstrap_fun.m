function [bs, ns, ms, logb] = bootstrap_fun (data, Nboot, n, m)

% call function, depending on whether or not m and n are specified
if nargin == 4
    bhat = bootstrp(Nboot, @(x) bootstrap_nmr_called(x, n, m), data);    % n,m can vary
elseif nargin == 3
    bhat = bootstrp(Nboot, @(x) bootstrap_nmr_called(x, n), data);    % m can vary
elseif nargin == 2
    bhat = bootstrp(Nboot, @(x) bootstrap_nmr_called(x), data); % m = 0, n = 2
else
    error('Not enough input arguments')
end

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