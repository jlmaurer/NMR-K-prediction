function [svs,h] = compute_semivariance(X,Z)
% X is the position vector and is n_dim x n_obs, Z is the vector of
% observation values and is n_obs x 1. 
[h, svs] = deal(zeros(length(Z))); 
for loop = 1:length(Z)
    for loop2 = 1:length(Z)
        svs(loop,loop2) = 0.5*(Z(loop) - Z(loop2))^2; 
        h(loop, loop2) = norm(X(loop)- X(loop2)); 
    end
end
svs = svs(:);
h = h(:); 
end