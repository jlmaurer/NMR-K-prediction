function [logKpred] = KGMfun(x, xdata)
% this function computes the likelihood for Kpredicted given K_measure

    [aa,bb] = size(x); 
    if aa ~=1 & bb~=1
        if bb > aa
            x = x'; 
        end
    end
    T2 = xdata(:,1);
    if size(xdata,2)==2
        phi = xdata(:,2);        
        lphi = log10(phi(:)); 
    end

    Temp = 20;  % temperature in degress C 
    rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
    eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
    Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds
    D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
    g = 9.8;    %m/s^2
    tort = 1/(1.5^2); 
    t1 = (rho(Temp)*g)/(8*eta(Temp)); % 
    
    num2 = 4*D(Temp)*Tb(Temp)*T2;
    denom2 = Tb(Temp) - T2; 

    if size(x,2) ==3
        ltinv = x(:,1); 
        B = x(:,2); 
        m = x(:,3); 
    elseif size(x,2) == 2
        ltinv = x(:,1); 
        B = x(:,2); 
        m = 0; 
    else
        ltinv = tort;
        B = x(1); 
        m = 0; 
    end
    f12 = (D(Temp)./B);  
    SQterm = sqrt(f12.^2 + (num2./denom2)); 

    % predicted data
    if size(xdata,2)==2
        logKpred = ltinv + log10(t1) + m*lphi + 2*log10(SQterm-f12);
    else
        logKpred = ltinv + log10(t1) + 2*log10(SQterm-f12);
    end

end