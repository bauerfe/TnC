%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% - only for biogeochemistry part - %%%%%%%%%%%%%%%%%%%%%
% - data from initial spinup must already be loaded in struct `prev_data` - 
%%%%%%%%%%%%% - data from bg spinup in struct `prev_bg` - %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: High vegetation ignored here!

% Run default PARAM file here, then overwrite plant- and nutrient related parameters with
% last values from previous spinup stage
run(initial_parameters_path);

%! Don't ignore bounds on nutrient reserves `rNc`, `rPc`, `rKc` in `Nutrients_Available`
OPT_IgnoreNutrientConcentrationBounds = 0;

% Mpar_L(1).NPK_res_ini = [35 5 15];  % [0 0 0];
Mpar_L(1).NPK_res_ini = [0 0 0];  % [0 0 0]; Don't provide initial NPK reserve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=prev_data.LAI_L(end,:);
B_L(1,:,:)= prev_data.B_L(end,:,:);
Rrootl_L(1,:)=prev_data.Rrootl_L(end,:);
PHE_S_L(1,:)=prev_data.PHE_S_L(end,:);
dflo_L(1,:)=prev_data.dflo_L(1,:);
AgeL_L(1,:)=prev_data.AgeL_L(end,:);
e_rel_L(1,:)=prev_data.e_rel_L(end,:);
hc_L(1,:) =prev_data.hc_L(end,:);
SAI_L(1,:)=prev_data.SAI_L(end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=prev_data.BLit(end);  %% [kg DM /m2 PFT] Litter Biomass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil Biogeochemistry part
Nreserve_L(1,:) = prev_data.Nreserve_L(end,:);
Preserve_L(1,:) = prev_data.Preserve_L(end,:);
Kreserve_L(1,:) = prev_data.Kreserve_L(end,:);

FNC_L(1,:) = 1;
NupI_L(1,:,:) = prev_data.NupI_L(end,:,:);
RexmyI(1,:,:) = prev_data.RexmyI(end,:,:);
NavlI(1,:,:) = prev_data.NavlI(end,:,:);

%% From spinup
P(1,:) = prev_bg;
P(P<0) = 0;

%%%%% Constant, BG-related parameters
Upl=0.01; %% Soil production [mm/yr]
HIST=0;
% Have to run `Biogeochemistry_IO` from root dir to ensure file path is correct
cd(root_dir);
[B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl);
B_IO.SC_par=[1 1 1 1];
PHs=7;

%%%%% Fertilizer Management
manure_data_path = [root_dir filesep 'data' filesep 'Rothamsted_manure_data.csv'];
B_IO = load_fertilizer(B_IO, dateNum1, NNd, crop_data, manure_data_path, OPT_Use_Fertilizer);

% Following values are mostly affeced by meteorological forcings, probably fine to stick with default values and not taking them from previous data
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=3;   % From Buckley et al., 2024 
Tdp(1,:)= 3*ones(1,ms);  % From Buckley et al., 2024 
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J �C/m^2 K]
TdpI_H(1,:)=3;  % From Buckley et al., 2024 
TdpI_L(1,:)=3;  % From Buckley et al., 2024 
%%% Snow_alb = soil_alb initial 
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
SP_wc(1)=0 ; %%[mm]
In_Litter(1)=0;
In_SWE(1)= 0;
ros(1)= 0;
t_sls(1)= 0;
e_sno(1) = 0.97;
tau_sno(1) = 0;
EK(1)=0; 
WAT(1) = 0;
ICE(1) = 0; 
IP_wc(1)=0; 
ICE_D(1)= 0;
FROCK(1)=0;
Ws_under(1)=1; 
%%%%%%%%%%%%%% Volume [mm]
O(1,:)=  Ofc;  
%%%%%%%%%%%%%%%%%%%%%%%%
%Vmin= 0*dz;
%Vmax= (Osat-Ohy)*dz;
%Vmax(ms) = Vmax(ms) + (Zbed-ZWT(1))*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%V(1,ms)= 0 + (Osat-Ohy)*(Zbed-ZWT(1));
%%%%%%