function [d, K, T2ML, phi, z, SumEch, log10K, log10T2, log10Porosity, SumEch_3s, SumEch_twm, SumEch_twm_3s] = loadnmrdata2(name)
    % this function loads a nmr data file namaed 'name' and defines the
    % variables. 
    
    if strcmp(name, 'A1') ==0 && strcmp(name, 'C1') == 0 && strcmp(name, 'gems_all') ==0 && strcmp(name, 'all_data') == 0
        in=(['..' filesep '..' filesep 'Data' filesep 'Aggregated_data' filesep name '.txt']);               %define the datafile name
        d=load(in);                         %load the raw data
        K= d(:,4);                %convert direct data to K in [m filesep s] 
        T2ML=d(:,2);                        %determine range of T2ML data, converts Vista Clara to [ms]
        phi=d(:,3);                         %determine range of NMR PHI data
        z=d(:,1);                           %depth
        SumEch = d(:,5);                    %Back out raw sum of echoes value from vista clara processed data
        SumEch_3s = d(:,6); 
        SumEch_twm = d(:,7); 
        SumEch_twm_3s = d(:,8); 
            
        log10K = log10(K); 
        log10T2 = log10(T2ML); 
        log10Porosity = log10(phi);
        
    elseif strcmp(name, 'A1') ==1 
        str = ['..' filesep 'Data' filesep 'dpnmr_site_' name];
        try
            load(str)
        catch
            keyboard;
        end
        T2ML = T2; 
        log10Porosity = log10(phi); 
        log10T2 = log10(T2ML);
        SumEch = 10.^SOE; 
        d = [z, T2ML, phi, K, SumEch];
        str2 = ['..' filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrA11.txt']; 
        try
            d2 = load(str2);
        catch
            keyboard; 
        end
        
        SumEch_3s1 = d2(:,6); 
        SumEch_twm1 = d2(:,7); 
        SumEch_twm_3s1 = d2(:,8); 
        str3 = ['..' filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrA12.txt'];
        d3 = load(str3);
        SumEch_3s = mean([SumEch_3s1, d3(:,6)], 2); 
        SumEch_twm = mean([SumEch_twm1, d3(:,7)], 2); 
        SumEch_twm_3s = mean([SumEch_twm_3s1, d3(:,8)], 2); 
        log10K = log10(K); 
        
    elseif strcmp(name, 'C1') == 1
        str = ['..' filesep 'Data' filesep 'dpnmr_site_' name]; 
        load(str)
        T2ML = T2; 
        log10Porosity = log10(phi); 
        log10T2 = log10(T2ML);
        SumEch = 10.^SOE; 
        d = [z, T2ML, phi, K, SumEch];
        str2 = ['..' filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrC1S.txt']; 
        d2 = load(str2);
        SumEch_3s1 = d2(:,6); 
        SumEch_twm1 = d2(:,7); 
        SumEch_twm_3s1 = d2(:,8); 
        str2 = ['..' filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrC1SW.txt']; 
        d2 = load(str2);
        SumEch_3s2 = d2(:,6); 
        SumEch_twm2 = d2(:,7); 
        SumEch_twm_3s2 = d2(:,8);
        str3 = ['..' filesep 'Data'  filesep 'Aggregated_data' filesep 'dpnmrC1SE.txt'];
        d3 = load(str3);
        SumEch_3s = mean([SumEch_3s1, SumEch_3s2, d3(:,6)], 2); 
        SumEch_twm = mean([SumEch_twm1, SumEch_twm2, d3(:,7)], 2); 
        SumEch_twm_3s = mean([SumEch_twm_3s1, SumEch_twm_3s2, d3(:,8)], 2); 
       
    elseif strcmp(name, 'gems_all') ==1
        
        str = ['..' filesep 'Data' filesep 'dpnmr_site_C1']; 
        load(str)
        str2 = ['..'  filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrC1S.txt']; 
        d2 = load(str2);
        SumEch_3s1 = d2(:,6); 
        SumEch_twm1 = d2(:,7); 
        SumEch_twm_3s1 = d2(:,8); 
        str2 = ['..' filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrC1SW.txt']; 
        d2 = load(str2);
        SumEch_3s2 = d2(:,6); 
        SumEch_twm2 = d2(:,7); 
        SumEch_twm_3s2 = d2(:,8);
        str3 = ['..' filesep 'Data'  filesep 'Aggregated_data' filesep 'dpnmrC1SE.txt'];
        d3 = load(str3);
        
        SumEch_3sC = mean([SumEch_3s1, SumEch_3s2, d3(:,6)], 2); 
        SumEch_twmC = mean([SumEch_twm1, SumEch_twm2, d3(:,7)], 2); 
        SumEch_twm_3sC = mean([SumEch_twm_3s1, SumEch_twm_3s2, d3(:,8)], 2); 
       
        zC = z; 
        T2C = T2; 
        Dkc = K; 
        phic = phi; 
        SumEchc = 10.^SOE; 

        clearvars -except zC T2C Dkc phic SumEchc SumEch_3sC SumEch_twmC SumEch_twm_3sC
        str = ['Data' filesep 'dpnmr_site_A1']; 
        load(str)
        T2ML = T2; 
        SumEch = 10.^SOE; 
        str2 = ['..' filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrA11.txt']; 
        d2 = load(str2);
        SumEch_3s1 = d2(:,6); 
        SumEch_twm1 = d2(:,7); 
        SumEch_twm_3s1 = d2(:,8); 
        str3 = ['..'  filesep 'Data' filesep 'Aggregated_data' filesep 'dpnmrA12.txt'];
        d3 = load(str3);
        SumEch_3sA = mean([SumEch_3s1, d3(:,6)], 2); 
        SumEch_twmA = mean([SumEch_twm1, d3(:,7)], 2); 
        SumEch_twm_3sA = mean([SumEch_twm_3s1, d3(:,8)], 2); 
       
        zA = z; 
        T2a = T2; 
        Dka = K; 
        phia = phi; 
        SumEcha = 10.^SOE;
        
        z = [zA; zC]; 
        K = [Dka; Dkc]; 
        T2ML = [T2a; T2C]; 
        SumEch = [SumEcha; SumEchc];
        phi = [phia; phic]; 
        SumEch_3s = [SumEch_3sA(:); SumEch_3sC(:)];
        SumEch_twm = [SumEch_twmA(:); SumEch_twmC(:)];
        SumEch_twm_3s = [SumEch_twm_3sA(:); SumEch_twm_3sC(:)];

        d = [z, T2ML, phi, K, SumEch, SumEch_3s, SumEch_twm, SumEch_twm_3s]; 
%         save([cd ' filesep bstrp_dat_VC_processed filesep gems_all.txt'], 'd', '-ascii')

        log10T2 = log10(T2ML); 
        log10Porosity = log10(phi); 
        log10K = log10(K); 

    elseif strcmp(name, 'all_data') == 1
        [~, Dkg, T2MLg, phig, zg, SumEchg, ~, ~, ~, SumEch_3sg, ...
            SumEch_twmg, SumEch_twm_3sg] = loadnmrdata2('gems_all');
        [~, Dklar, T2MLlar, philar, zlar, SumEchlar, ~, ~, ~, SumEch_3slar, ...
            SumEch_twmlar, SumEch_twm_3slar] = loadnmrdata2('dpnmr_larned_all');
        [~, Dkleq, T2MLleq, phileq, zleq, SumEchleq, ~, ~, ~, SumEch_3sleq, ...
            SumEch_twmleq, SumEch_twm_3sleq] = loadnmrdata2('dpnmr_leque_all');
        

        z = [ zg; zlar; zleq]; 
        K = [Dkg; Dklar; Dkleq]; 
        T2ML = [T2MLg; T2MLlar; T2MLleq]; 
        SumEch = [SumEchg; SumEchlar; SumEchleq];
        phi = [phig; philar; phileq];
        SumEch_3s = [SumEch_3sg; SumEch_3slar; SumEch_3sleq];
        SumEch_twm = [SumEch_twmg; SumEch_twmlar; SumEch_twmleq];
        SumEch_twm_3s = [SumEch_twm_3sg; SumEch_twm_3slar; SumEch_twm_3sleq];
        d = [z, T2ML, phi, K, SumEch, SumEch_3s, SumEch_twm, SumEch_twm_3s]; 
        
        log10T2 = log10(T2ML); 
        log10Porosity = log10(phi); 
        log10K = log10(K); 
        
    else
        error('No such data file')
    end 
end