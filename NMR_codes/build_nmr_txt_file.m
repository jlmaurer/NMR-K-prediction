function [z, T2ML, phi, Dk, soe, soe_3s, soe_twm, soe_twm_3s,k_SDR_VC, t, S, CapH20, FreeH20, k_TC] = build_nmr_txt_file(name)
% convert raw data files to text files in format to be used by
% NMR_script_main

% read in data files
in1 = [pwd '\bstrp_dat_VC_raw\' name '\' name '_SE_decay' '.txt']; 
in2 = [pwd '\bstrp_dat_VC_raw\' name '\' name '.txt'];
in3 = [pwd '\bstrp_dat_original\' name '.txt'];

decaycurve = load(in1); 
dparam = load(in2); 
olddat = load(in3); 

% set needed variables
z = dparam(:,1); 

t = decaycurve(1,:); 
S = decaycurve(2:end, :); 

Dk = olddat(:,4)*1.16e-5; 
z_dk = olddat(:,1); 

%% Match depths of permeability measurements from my old data
% with the depths in the new processed data -- need to find nearest
% neighbors, exlude NMR points that do not have corresponding k
% measurement. 

for i = 1:length(Dk)
    [errorz(i), ind(i)] = min(abs(z_dk(i) - z)); 
end

%% Set rest of variables
phi = dparam(ind,2); 
CapH20 = dparam(ind,3); 
FreeH20 = dparam(ind,4); 
T2ML = dparam(ind,5); 
k_SDR_VC = dparam(ind,6); 
k_TC = dparam(ind,7); 
soe = dparam(ind,8); 
soe_3s = dparam(ind,9); 
soe_twm = dparam(ind,10); 
soe_twm_3s = dparam(ind,11); 

z = z(ind); 

%% check SOE and T2ML values from old data and new processed data. 


%%%%%%%%%%% Plot actual decay curves at each depth
    figure;
    plot(t, S)
    xlabel('time')
    ylabel('S(t)')
    legend(num2str(z))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5    

% create text file output for use with NMR_script_main
datamat = [z(:), T2ML(:), phi(:), Dk(:), ...
    soe, soe_3s(:), soe_twm(:), soe_twm_3s(:)]; 

save([cd '/bstrp_dat_VC_processed/' name '.txt'], 'datamat', '-ascii')
end