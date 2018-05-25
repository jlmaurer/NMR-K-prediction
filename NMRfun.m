function [loglike, kpred, likelihood] = NMRfun(x, Dk, phi, T2, ip)
% this function computes the likelihood for Kpredicted given K_measure
    if length(x) ==5
        T2B = x(1); 
        logb = x(2); 
        n = x(3); 
        m = x(4); 
        sig = x(5); 
    elseif length(x) == 3
        logb = x(1); 
        n = x(2); 
        sig = x(3); 
    else
        error('Bad x')
    end
    
    % put in log space
    kk = log10(Dk); 
    lphi = log10(phi); 
    
    % compute forward model for chosen model type
    if ip == 1
        kpred = logb + m*lphi - n*log10((T2.^-1)-(T2B.^-1));    % with porosity
    elseif ip == 2
        kpred = logb - n*log10((T2.^-1)-(T2B.^-1));    % without porosity
    elseif ip == 0
        kpred = logb + n*log10(T2);                    % without T2B
    elseif ip == 3
        kpred = logb*(T2.^n); 
    end
    
    r = kk - kpred;
    tau = (1./sig.^2); 
    n = length(r); 

   liketerm1 = (-.5*tau)*(r'*r); 
   liketerm2 = (n/2)*log(tau) - (n/2)*log(2*pi);
   loglike = liketerm2 + liketerm1; 
   likelihood = exp(loglike);
end
