%! The actual spinup
function [B, R_litter, R_litter_sur, R_microbe, R_bacteria, R_ew, VOL, BfixN, Min_N, Min_P, RmycAM, RmycEM, ...
    N2flx, NH4_Uptake, NO3_Uptake, P_Uptake, K_Uptake, LEAK_NH4, LEAK_NO3, LEAK_P, LEAK_K, LEAK_DOC, ...
    LEAK_DON, LEAK_DOP, Lk] = biogeochemistry_spinup(Lat, Lon, N_epochs, Nyears_epoch, ISOIL, Zbio, rsd, ...
    Pcla, Psan, PH, Zs, Ts, Ta, Psi_s, Se, Se_fc, V, VT, T_L, T_H, Lk, RexmyI, Ccrown, B, rtol, atol, save_path, ...
    fertilizer_ref_date, crop_data, manure_data_path, OPT_Use_Fertilizer)
    
    ZBIOG = 0.001 * Zbio; %%[m] %! Biogeochemically active layer depth in m?
    rsd = mean(rsd); %% density dry soil [kg/m^3]
    Zs = 0.001 * max(Zs); %% [m] max. soil depth?
    
    Ta = nanmean(Ta);
    T = nanmean(T_L + T_H) * 24;
    Lk = mean(Lk) * 24;
    RexmyI = mean(RexmyI, 1); %% [gC m2 / day] %% mean(RexmyI,1)
    Broot = mean(B(:, 1, 3));

    IS = squeeze(mean(ISOIL, 1));

    %! These initial conditions are taken from another file and are assumed to be arbitrary
    %%%%%%% CARBON POOL %%%%%%%%%%%
    Btm1(1)=121; %%% B1 Above-ground Litter Metabolic
    Btm1(2)=315; %%% B2 Above-ground Litter Structural - Cellulose/Hemicellulose
    Btm1(3)=51.5; %%% B3 Above-ground Litter Structura - Lignin
    Btm1(4)=0;  % 3015; %%% B4 Above-ground Woody  - Cellulose/Hemicellulose
    Btm1(5)=0; %%% B5 Above-ground Woody - Lignin
    Btm1(6)=47.1; %%% B6 Below-ground Litter Metabolic
    Btm1(7)=381.8; %%% B7 Below-ground Litter Structural - Cellulose/Hemicellulose
    Btm1(8)=95.0; %%% B8 Below-ground Litter Structura - Lignin
    Btm1(9)= 500;  % 3077; %%% B9  SOM-POC- lignin
    Btm1(10) = 1144; %%% B10 SOM-POC -Cellulose/Hemicellulose
    Btm1(11) = 8631; %%% B11 SOM-MOC 
    Btm1(12) = 10.21; %%%B12 DOC - for bacteria 
    Btm1(13) = 8.17; %B13 DOC - for fungi 
    Btm1(14) = 0.069; %%% B14 Enzyme for decomposition of POC-Bact 
    Btm1(15) = 0.046; %%% B15 Enzyme for decomposition of POC-Fung
    Btm1(16) = 0.028; %% B16 Enzyme for decomposition of MOC-Bact
    Btm1(17) = 0.077; %%% B17 Enzyme for decomposition of MOC-Fung
    Btm1(18) = 32.96; %%% B18 Bacteria pool
    Btm1(19) = 126.0; %%%  B19 Fungi saprotrophic 
    Btm1(20) = 33.35; %%% B20 AM-Mycorrhizal - C 
    Btm1(21) = 0; %%%  B21 EM-Mycorrhizal - C 
    Btm1(22) = 1.00; %%%  B22 Earthworms - C 
    %%
    %%%%%% NITROGEN POOL
    Btm1(23) = 4.83; %%% B23 Nitrogen Above-ground Litter
    Btm1(24) = 0; % 20.8;%%% B24 Nitrogen Above-ground Woody
    Btm1(25) = 4.13; %%% B25 Nitrogen Below-ground Litter
    Btm1(26) = 936.7; %%% B26 Nitrogen SOM
    Btm1(27) = 6.35;%%% B27 Nitrogen Bacteria 
    Btm1(28) = 19.40;%%% B28 Nitrogen Fungi 
    Btm1(29) = 1.85;%%% B29 AM Mycorrhizal - N 
    Btm1(30) = 0;%%% B30 EM Mycorrhizal - N  
    Btm1(31) = 0.282;%%% B31 Nitrogen Ione Ammonium NH4+
    Btm1(32) = 0.029;%%%B32 Nitrogen Nitrate NO3-
    Btm1(33) = 0.013; %%% B33 DON
    Btm1(34) = 0.098; %%% B34 Earthworms - N 
    %%%%%% PHOSPHORUS POOL
    Btm1(35) = 0.345; %%% B35 phosphorus Above-ground Litter
    Btm1(36) = 0; % 1.487; %%% B36 phosphorus Above-ground Woody
    Btm1(37) = 0.295; %%% B37 phosphorus Below-ground Litter
    Btm1(38) = 168.8; %%% B38 phosphorus SOM
    Btm1(39) = 2.05;%%%%%% B39 phosphorus Bacteria 
    Btm1(40) = 3.15;%%%% B40 phosphorus Fungi 
    Btm1(41) = 0.277;%%% B41 AM - Mycorrhizal - P
    Btm1(42) = 0.0;%%%% B42 EM - Mycorrhizal - P
    Btm1(43) = 0.0199;%%% B43 phosphorus Mineral
    Btm1(44) = 150;%%%%% B44 phosphorus primary
    Btm1(45) = 0.397;%% B45 phosphorus secondary
    Btm1(46) = 15;%%%%% B46 phosphorus occluded
    Btm1(47) = 0.0001; %%%% B47 DOP
    %%%%%% POTASSIUM POOL
    Btm1(48) = 2.41;%%%% B48 Potassium Above-ground Litter
    Btm1(49) = 0; % 10.41; %%% B36 Potassium  Above-ground Woody
    Btm1(50) = 0.637;%%%% B50 Potassium  Below-ground Litter
    Btm1(51) = 11.77; %%%% B51 Potassium SOM
    Btm1(52) = 0.1168; %%%% B52 Potassium  Mineral  solution
    Btm1(53) = 0.103; %%%% B53 Potassium  exchangeable
    Btm1(54) = 2.332; %%% B54 Potassium fixed or non-exchangeable 10- 80
    Btm1(55) = 502.9; %%%% B55  Potassium in the lattice of certain primary minerals  130-9500
    %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% External Inputs
    %DepN=  0.003287; %   gN/m2 day
    %DepP= 2.8385e-006; %  gP/m2 day
    %DepK= 9.698e-005; %%  gK/m2 day
    %FertN=0*ones(1,366);
    %FertP=0*ones(1,366);
    %FertK=0*ones(1,366);
    Upl=0.01; %% Soil production [mm/yr]
    HIST=0; 
    %[B_IO]=Biogeochemistry_IO(Zs*1000,Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs*1000,Lat,Lon,Upl);
    ExEM = 0; 
    B_IO.SC_par= [1 1 1 1]; %%% (Maint Bact , Maint Fungi,  EP, Macrofuanal)

    if OPT_Use_Fertilizer
        disp("Including fertilizer in BG-spinup")
        B_IO = load_fertilizer(B_IO, fertilizer_ref_date, 365, crop_data, manure_data_path, true);
    end

    FertN=mean(B_IO.FertN);
    DepN=B_IO.DepN;
    FertP=mean(B_IO.FertP);
    DepP=B_IO.DepP; 
    FertK=mean(B_IO.FertK); 
    DepK=B_IO.DepK; 
    Tup_P=B_IO.Tup_P; 
    Tup_K=B_IO.Tup_K; 
    SC_par = B_IO.SC_par; 

    % Manure
    N_Man =B_IO.N_Man; %% Manure  [gC/gN]
    P_Man =B_IO.P_Man; % Manure  [gC/gP]
    K_Man =B_IO.K_Man;% Manure  [gC/gK]

    Lig_fr_Man =B_IO.Lig_fr_Man; %% Lignin fraction in Manure  [g Lignin / g DM] 
    frac_to_metabolic_Man = 0.85 - 0.018*(N_Man*2*Lig_fr_Man); 
    frac_to_metabolic_Man(frac_to_metabolic_Man<0)=0;
    ManF = mean(B_IO.ManF);  %%% [gC /m2 day] 
    IMAN(1)= frac_to_metabolic_Man*(ManF) ; %% met_sur_lit  [gC/m^2 d]
    IMAN(2)= (1-frac_to_metabolic_Man)*(ManF)*Lig_fr_Man ;%%  str_sur_lit_lig [gC/m^2 d]
    IMAN(3)=  (1-frac_to_metabolic_Man)*(ManF)*(1-Lig_fr_Man);  %% str_sur_lit_nlig [gC/m^2 d]
    IMAN(4)=  ManF./N_Man; % [gN/m^2 day]  
    IMAN(5)=  ManF./P_Man; % [gP/m^2 day]
    IMAN(6)=  ManF./K_Man; % [gK/m^2 day]

    opt_cons_CUE=1;
    [BiogeoPar]=Biogeochemistry_Parameter(opt_cons_CUE);
    
    Ndays_epoch = 365 * Nyears_epoch;
    t=1:Ndays_epoch;

    %! We don't need to reset these values between epochs, as the values in each interation are independent of previous values
    B=zeros(Ndays_epoch,55);
    R_litter=zeros(Ndays_epoch,1);
    R_litter_sur=zeros(Ndays_epoch,1);
    R_microbe=zeros(Ndays_epoch,1);
    R_ew=zeros(Ndays_epoch,1);
    VOL=zeros(Ndays_epoch,1);
    BfixN = zeros(Ndays_epoch,1);
    Min_N = zeros(Ndays_epoch,1);
    Min_P = zeros(Ndays_epoch,1);
    R_bacteria= zeros(Ndays_epoch,1);
    RmycAM = zeros(Ndays_epoch,1);
    RmycEM = zeros(Ndays_epoch,1);
    N2flx = zeros(Ndays_epoch,1);
    NH4_Uptake=zeros(Ndays_epoch,1);
    NO3_Uptake = zeros(Ndays_epoch,1);
    P_Uptake=zeros(Ndays_epoch,1);
    K_Uptake = zeros(Ndays_epoch,1);
    LEAK_NH4 = zeros(Ndays_epoch,1);
    LEAK_NO3 = zeros(Ndays_epoch,1);
    LEAK_P = zeros(Ndays_epoch,1);
    LEAK_K = zeros(Ndays_epoch,1);
    LEAK_DOC = zeros(Ndays_epoch,1);
    LEAK_DON = zeros(Ndays_epoch,1);
    LEAK_DOP = zeros(Ndays_epoch,1);
    
    for nn=1:N_epochs
        for j=1:Ndays_epoch
            %%%%%%%%%
            [LEAK_NH4(j),LEAK_NO3(j),LEAK_P(j),LEAK_K(j),LEAK_DOC(j),LEAK_DON(j),LEAK_DOP(j)]= Biogeo_Leakage(Btm1,Lk,V,BiogeoPar);
            %%%%
            [NH4_Uptake(j),NO3_Uptake(j),P_Uptake(j),K_Uptake(j)]= Biogeo_uptake(Btm1,Broot,Ts,T,VT,Ccrown,ExEM,BiogeoPar);
            %NH4_Uptake(j)=0; NO3_Uptake(j)=0; P_Uptake(j)=0; K_Uptake(j)=0;
            %%%%
            [BfixN(j)]= Biogeo_Bio_fixation(0,0,Btm1,RexmyI,Ts);
            %
            %%%

            [dB,R_litter(j),R_microbe(j),R_litter_sur(j),R_ew(j),VOL(j),N2flx(j),Min_N(j),Min_P(j),R_bacteria(j),RmycAM(j),RmycEM(j)]= BIOGEOCHEMISTRY_DYNAMIC3(t(j),Btm1,ZBIOG,rsd,IS,Ts,Ta,Psi_s,PH,Se,Se_fc,FertN,DepN,BfixN(j),FertP,DepP,FertK,DepK,...
                NH4_Uptake(j),NO3_Uptake(j),P_Uptake(j),K_Uptake(j),LEAK_DOC(j),LEAK_NH4(j),LEAK_NO3(j),LEAK_P(j),LEAK_K(j),LEAK_DON(j),LEAK_DOP(j),Tup_P,Tup_K,ExEM,Pcla,Psan,BiogeoPar,SC_par,IMAN,opt_cons_CUE); 
        
            %%%%
            if isreal(sum(dB))==0 || isnan(sum(dB)) == 1
                disp(['Stopping in epoch ' nn ', day ' jj ' because of some break condition related to variable `dB`'])
                break
            end
            %%%
            %%%%
            B(j,:)=Btm1+dB;
            Btm1=B(j,:);
            %%%
        end 

        if ~strcmp(save_path, 'None')
            save([save_path filesep sprintf('%04d', nn) '.mat']);
        end

        % Compare how many values of B have changed by a certain relative amount during the epoch
        B(B<1e-10) = 1e-10;
        diff_B = abs((B(end,:) - B(1, :)));
        rel_diff_B = diff_B ./ B(end, :);
        n_ss = sum((rel_diff_B < rtol) | (diff_B < atol));
        disp([ 'Spin-up epoch ' num2str(nn) ': ' num2str(n_ss) ' of 55 elements in steady state.'])

        if n_ss == 55
            disp("Spin-up complete")
            break
        end
    end

    %! To load data across full spinup do the following
    % B_full = zeros(Ndays_epoch * nn, 55);
    % for i=1:nn
    %     bg_curr=load([save_path filesep sprintf('%04d', i) '.mat']);
    %     B_full(1+(i-1)*Ndays_epoch : i*Ndays_epoch, :) = bg_curr.B;
    % end
    
    
    
    
    %%%%%%%% SOIL BIOGEOCHEMISTRY BALANCE CHECK
    IS=IS*Ndays_epoch;
    C_exp = sum(IS(1:9));
    N_exp = sum(IS(10:12)) + DepN*Ndays_epoch + FertN*Ndays_epoch;
    P_exp = sum(IS(13:15)) + DepP*Ndays_epoch + Tup_P*Ndays_epoch + FertP*Ndays_epoch ;
    K_exp = sum(IS(16:18)) + DepK*Ndays_epoch + Tup_K*Ndays_epoch + FertK*Ndays_epoch ;
    dP_soil = B(1,:) - B(end,:);
    C_out = sum(LEAK_DOC)+sum(R_litter)+sum(R_microbe)+sum(R_ew);
    N_out = sum(sum(NH4_Uptake +NO3_Uptake)) + sum(LEAK_NH4) + sum(LEAK_NO3) + sum(LEAK_DON) + sum(VOL) +sum(N2flx);
    P_out = sum(sum(P_Uptake)) + sum(LEAK_DOP) + sum(LEAK_P);
    K_out = sum(sum(K_Uptake)) +sum(LEAK_K);
    
    
    CkC_s = sum(dP_soil([1:22])) + C_exp - C_out;
    CkN_s= sum(dP_soil([23:34]))+ N_exp - N_out;
    CkP_s= sum(dP_soil([35:47]))+ P_exp - P_out;
    CkK_s= sum(dP_soil([48:55]))+ K_exp - K_out;
end
