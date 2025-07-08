%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [wheat, barley,grass, Potatoes, RS, Peas]

% MAIN PARAMATERS TO CHANGE
% I have played around with the below values,largely out of curiosity, double check 
% values are realistic... 

Sl_L = [0.03 0.032 0.04 0.024 0.03 0.029]; % 
Nl_L= [15 21 30 25 18 22]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
age_cr_L= [95 90 110 100 85 58]; %% [day] Critical Leaf Age
dmg_L= [70 60 65 70 45 30]; %%% Tree 30 Grasses Day of Max Growth

%In order to play with yields, the parameter soCrop_L, and ff_r_L are good targets


Vmax_L = [120 120 70 75 118 50]; %  % Wheat, barley, % [160 150 125 80 119 58] before paying with yields...
%Usually, Vmax < 150, even though crops have the highest values and I also see that we need very high values for certain crops. 


%Realistic values are Sl < 0.055, crop have high values, but anything above 0.050 is extreme. 

%Realistic values are Sl < 0.055, crop have high values, but anything above 0.050 is extreme. 
r_L = [0.016 0.025 0.021 0.025 0.025 0.03];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ] % was 0.25
gR_L= [0.16 0.22 0.21 0.25 0.25 0.3]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%age_cr_L= [100 92 120 100 120 78];%[120 125 250 90 120 90]; %% [day] Critical Leaf Age
Tlo_L = [5.5 4 5 7 6.7 7.5]; %% Mean Temperature for Leaf onset
%dmg_L= [71 50 90 80 80 35]; %%% Tree 30 Grasses Day of Max Growth
mjDay_L = [250 250 300 250 200 250]; %% Maximum Julian day for leaf onset % so they dont come up in december... 
LDay_min_L =[10 10.5 4 10 7.2 9]; %% Minimum Day duration for leaf onset % 8 worked
LDay_cr_L = [13.0 11 12 12.4 12 12.0]; %%%  Threshold for senescence day light [h] [13.0 12.1 10 12.4 12 12.0];
Sl_emecrop_L= [0.015 0.018 0.02 0.015 0.025 0.005]; %%% Additional SLA at emergence 
MHcrop_L =[1.2 1.2 0.7 0.8 1.2 0.7]; %%[m] maximum crop height 
ff_r_L= [0.5 0.5 0.5 0.50 0.50 0.6]; % fraction allocated to fruit ? 

% one source states that 80% of field pea yield is off the main stem. 

DSE_L =[0.649 0.649 0.649 0.649 0.549 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72 72 72 72 72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent % at high temp, increase for more photosynthesis, oposite for at low temps. 

%Mpar_L(6).Crop_B=[8 11];  % second one is C in seed Remember we can play
%with these too ! 

%
Pcla= 0.40;
Psan= 0.10;
Porg= 0.025;
Psi_sto_50_L = [-2.9 -2.9  -2.9  -2.9  -2.9  -2.9 ];%  %% [MPa]  Water Potential at 50% loss conductivity
Trr_L = [6.5 6.5 6.5 6.5 6.5 6.5]; %% Translocation rate [gC /m^2 d]

ZR95_L = [950 900 600 500 600 500]; %% [mm] Depth of root
gcoef_L = [3.6 3.6 3.6 3.6 3.6 3.6]; % [gC/m2 day]
LtR_L = [1.0 1.1 1.2 1.0 1.2 1.6]; %%% Leaf to Root ratio maximum
soCrop_L = [0.10 0.10 0.1 0.10 0.10 0.10]; % no idea what this does !!!!! should affect yield

%%
% affects sensitivity of Vcmax to temperature. 

%16/10/2003 -> 04/08/2004 -> Winter Wheat [1]
%29/9/2004 -> 14/7/2005 -> Winter Barley [2]
%09/09/2005 -> Senesence -> Grass Cover crops [3] - 
%05/05/2006 -> Potatoes -> Dammaged Crop no harvest [4]
%19/10/2006-> 15/7/2007 -> Winter Wheat [5] - no rain...added
%28/8/2007 -> 16/7/2008 -> Rape Seed [6]
%7/10/08 -> 21/07/09 -> Winter Wheat [7] - no work
%12/08/09 -> Senesence -> Cover crop [8] 
%9/05/10 -> 19/7/10 -> Peas [9]
%15/10/10 -> 2/8/11 -> Winter Wheat [10]
%4/9/11 -> 9/7/2012 -> Winter Barley [11] %%%%%%%%%
% 4/9/13 -> 28/7/13 -> Rape Seed [12]
%19/10/13 ->  24/7/14 -> Winter Wheat [13] - AFter this are the new ones
%29/09/14- 4/7/15 -> Winter Barley [14]
%3/8/15 -> sensence cover crop [15]
%09/05/16 -> 25/7/16  peas [16]
%12/10/16 -> 19/07/17 - Winter Wheat [17]
%30/08/17 -> 12/07/18 -Rape Seed [18]
%11/10/18 -> 19/07/19 -> Winter Wheat [19] - 6 new crops in total 

%wheat, barley,grass, Potatoes, RS, Peas.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%
cur_dir=cd;
cd(Directory)
%%%%%%%%%%%%%%%%
%%% Rainfall Disaggregation
a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%
fpr=1;
SvF=1; %% Sky View Factor
SN=0; %% Stream Identifier
Slo_top=0;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared = 1; %% 1-Frock
aR =1; %%% anisotropy ratio
%Kh=Ks*aR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot = NaN; %% [mm/h] Conductivity at the bedrock layer
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm = 4.0; %% Reference Height (3.0 m when the crop is less than 1 m)
%%%%%%%%%%%%%%%%%%
%%%%%%% VEG. SPECIES  --- Low Grasses -- High Decidous
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare =  1.0; Ccrown = [0.0 0.0 0.0 0.0 0.0 0.0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL INPUT %%%%  silty clay loam

Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%%%%%%%%%%%%%%
clear Oss_Lp Owp_Lp
%%%%%%%%%%%%%%
Oice = 0;
%%%%%%%%%%%%%
Zs = [0 100 200 500 800 1000 1250 1500 2500 4000 5000]; %%% [ms+1]
Zs = [0 10 20 40 60 1000 1250 1500 2500 4000 5000]; %%% [ms+1] % to make it compatable with OBS SWC
Zs = [0 100 200 500 800 1000 1250 1500 2500 4000 5000]; %%% [ms+1]
Zs = [0 10 20 50 100 250 500 1000 1500 2500 5000]; % NOW USE THIS FOR ALL SITE RUNS


Zdes = 10;
Zinf=  10;
Zbio = 250;
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); %%% Infiltration Depth Layer fraction
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms);
for i = 1:ms
    if i>1
        Dz(i)= (dz(i)+ dz(i-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(i)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5;
In_max_rock=0.1; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=-0.8;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
%%%%%
ZR95_H = [0 0 0 0 0 0]; %% [mm]
%ZR95_L = [0 1200]; %% [mm]
%ZR95_L = [900 1200]; %% [mm]
ZR50_H = [NaN NaN NaN NaN NaN NaN];
ZR50_L = [NaN NaN NaN NaN NaN NaN];
ZRmax_H = [NaN NaN NaN NaN NaN NaN];
ZRmax_L = [NaN NaN NaN NaN NaN NaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2 0.2 0.2 0.2 0.2 0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.2 0.2 0.2 0.2 0.2 0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [NaN NaN NaN NaN NaN NaN]; %%[cm]
d_leaf_L= [5.0 5.0 5.0 5.0 5.0 5.0];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=[0.2 0.2 0.2 0.2 0.2 0.2]; %%% Canopy Nitrogen Decay
KnitL=[0.15 0.15 0.15 0.15 0.15 0.15] ; %%% Canopy Nitrogen Decay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mSl_H = [0.0 0.0 0.0 0.0 0.0 0.0];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [0.0 0.0 0.0 0.0 0.0 0.0]; % 0.001; %% [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[0.081 0.081 0.081 0.081 0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000 1000 1000 1000 1000 1000]; %%[Pa]
a1_H=[7 7 7 7 7 7];
go_H=[0.01 0.01 0.01 0.01 0.01 0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[4 4 4 4 4 4]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649 0.649 0.649 0.649 0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72 72 72 72 72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf Inf Inf Inf Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[1.97 1.97 1.97 1.97 1.97 1.97]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.081 0.081 0.081 0.081 0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons] % note that f ϵ = 0.081 [µmolCO2 µmol−1 photons]
%for C3 and ϵ = 0.040 [µmolCO2 µmol−1 photons] for C4 plants are typically used

Do_L=[1000 1000 1000 1000 1000 1000]; %%[Pa]
a1_L=[8 8 8 8 8 8];
go_L=[0.01 0.01 0.01 0.01 0.01 0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3 3 3 3 3 3 ];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants

gmes_L=[Inf Inf Inf Inf Inf Inf];
rjv_L= [2.4 2.4 2.4 2.4 2.4 2.4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pss_H = [800]; Pwp_H = [3000]; %%% [kPa]
%Pss_L = [500]; Pwp_L = [3500]; %%% [kPa]
Psi_sto_50_H =  [-2.0 -2.0 -2.0 -2.0 -2.0 -2.0] ;%% [MPa]  Water Potential at 50% loss conductivity
Psi_sto_00_H =  [-0.5 -0.5 -0.5 -0.5 -0.5 -0.5]; %% [MPa]  Water Potential at 2% loss conductivity
%%% Leaf
PsiL00_H =  [-2.5 -2.5 -2.5 -2.5 -2.5 -2.5] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  [-3.5 -3.5 -3.5 -3.5 -3.5 -3.5]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = [ 5 5 5 5 5 5] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200 1200 1200 1200 1200 1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [15.0 15.0 15.0 15.0 15.0 15.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000 80000 80000 80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-3.5 -3.5 -3.5 -3.5 -3.5 -3.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150 150 150 150 150 150]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%
Psi_sto_00_L = [-0.7 -0.7 -0.7 -0.7 -0.7 -0.7];%  %% [MPa]  Water Potential at PLCs loss conductivity
%%% Leaf
PsiL00_L =  [-1.2 -1.2 -1.2 -1.2 -1.2 -1.2]; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L =  [-3.5 -3.5 -3.5 -3.5 -3.5 -3.5];%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [5 5 5 5 5 5 ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200 1200 1200 1200 1200];  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [0.0 0.0 0.0 0.0 0.0 0.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000 80000 80000 80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-4.5 -4.5 -4.5 -4.5 -4.5 -4.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150 150 150 150 150 150]; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 

%%%% Growth Parameters
PsiG50_H= [-0.45 -0.45 -0.45 -0.45 -0.45 -0.45];  %%[MPa]
PsiG99_H= [-1.2 -1.2 -1.2 -1.2 -1.2 -1.2];  %%[MPa]
gcoef_H = [3.5 3.5 3.5 3.5 3.5 3.5]; % [gC/m2 day]
%%------
PsiG50_L= [-1.2 -1.2 -1.2 -1.2 -1.2 -1.2];
PsiG99_L= [-3.5 -3.5 -3.5 -3.5 -3.5 -3.5];

%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(0);
[PFT_opt_L(1)]=Veg_Optical_Parameter(16);
[PFT_opt_H(2)]=Veg_Optical_Parameter(0);
[PFT_opt_L(2)]=Veg_Optical_Parameter(16);
[PFT_opt_H(3)]=Veg_Optical_Parameter(0);
[PFT_opt_L(3)]=Veg_Optical_Parameter(16);
[PFT_opt_H(4)]=Veg_Optical_Parameter(0);
[PFT_opt_L(4)]=Veg_Optical_Parameter(16);
[PFT_opt_H(5)]=Veg_Optical_Parameter(0);
[PFT_opt_L(5)]=Veg_Optical_Parameter(16);
[PFT_opt_H(6)]=Veg_Optical_Parameter(0);
[PFT_opt_L(6)]=Veg_Optical_Parameter(16);

OM_H=[1 1 1 1 1 1];
OM_L=[1 1 1 1 1 1];
%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.016 0.016 0.016 0.016 0.016 0.016]; % [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [30 30 30 30 30 30]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
[Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
[Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
[Stoich_H(4)]=Veg_Stoichiometric_Parameter(Nl_H(4));
[Stoich_H(5)]=Veg_Stoichiometric_Parameter(Nl_H(5));
[Stoich_H(6)]=Veg_Stoichiometric_Parameter(Nl_H(6));
%PLNR_H = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_H = [0.030 0.030 0.030 0.030 0.030 0.030];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25 0.25 0.25 0.25 0.25 0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [1 1 1 1 1 1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [1/365 1/365 1/365 1/365 1/365 1/365]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [2/365 2/365 2/365 2/365 2/365 2/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [7 7 7 7 7 7]; %% [°C] Cold Leaf Shed
drn_H=  [1/1095 1/1095 1/1095 1/1095 1/1095 1/1095]; %% turnover root  [1/d]
dsn_H= [1/365 1/365 1/365 1/365 1/365 1/365]; % normal transfer rate sapwood [1/d]
age_cr_H= [150 150 150 150 150 150]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.95 0.95 0.95 0.95 0.95 0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN NaN NaN NaN NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [ 12.9 12.9 12.9 12.9 12.9 12.9]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN NaN NaN NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN NaN NaN NaN NaN NaN]; 
dmg_H= [35 35 35 35 35 35]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.01 0.01 0.01 0.01 0.01 0.01];
Trr_H = [3.5 3.5 3.5 3.5 3.5 3.5]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180 180 180 180 180 180]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.58 12.58 12.58 12.58 12.58 12.58]; %% Minimum Day duration for leaf onset
LtR_H = [0.8 0.8 0.8 0.8 0.8 0.8]; %%% Leaf to Root ratio maximum
Mf_H= [1/50 1/50 1/50 1/50 1/50 1/50]; %% fruit maturation turnover [1/d]
Wm_H= [0 0 0 0 0 0 ] ; % wood turnover coefficient [1/d]
eps_ac_H = [1 1 1 1 1 1]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [12.30 12.3 12.30 12.3 12.30 12.3]; %%%  Threshold for senescence day light [h]
Klf_H =[1/15 1/15 1/15 1/15 1/15 1/15]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74 0.74 0.74 0.74 0.74 0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26 0.26 0.26 0.26 0.26 0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1 0.1 0.1 0.1 0.1 0.1]; %% Reference allocation to Fruit and reproduction
soCrop_H = [NaN NaN NaN NaN NaN NaN]; 
MHcrop_H =[NaN NaN NaN NaN NaN NaN]; %% Maximum height crop  
Sl_emecrop_H = [NaN NaN NaN NaN NaN NaN]; %%% Additional SLA at emergence 
[ParEx_H(1)]=Exudation_Parameter(0); 
[ParEx_H(2)]=Exudation_Parameter(0); 
[ParEx_H(3)]=Exudation_Parameter(0);
[ParEx_H(4)]=Exudation_Parameter(0);
[ParEx_H(5)]=Exudation_Parameter(0);
[ParEx_H(6)]=Exudation_Parameter(0);
[Mpar_H(1)]=Vegetation_Management_Parameter;
[Mpar_H(2)]=Vegetation_Management_Parameter;
[Mpar_H(3)]=Vegetation_Management_Parameter;
[Mpar_H(4)]=Vegetation_Management_Parameter;
[Mpar_H(5)]=Vegetation_Management_Parameter;
[Mpar_H(6)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nl_L= [16 16 30 20 16 17]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
[Stoich_L(2)]=Veg_Stoichiometric_Parameter(Nl_L(2));
[Stoich_L(3)]=Veg_Stoichiometric_Parameter(Nl_L(3));
[Stoich_L(4)]=Veg_Stoichiometric_Parameter(Nl_L(4));
[Stoich_L(5)]=Veg_Stoichiometric_Parameter(Nl_L(5));
[Stoich_L(6)]=Veg_Stoichiometric_Parameter(Nl_L(6));
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]

%LAI_max_L= [2.0 10.0];2
aSE_L= [5 5 5 5 5 5]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -
dd_max_L= [1/50 1/50 1/50 1/50 1/50 1/50];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L = [7/365 7/365 7/365 7/365 7/365 7/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [0.0 0.0 0.0 0.0 0.0 0.0]; %% [°C] Cold Leaf SLed
drn_L=  [1/365 1/365 1/365 1/365 1/365 1/365]; %% turnover root  [1/d]
dsn_L= [0 0 0 0 0 0]; % normal transfer rate sapwood [1/d]
Bfac_lo_L= [0.99 0.99 0.99 0.99 0.99 0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN NaN NaN NaN NaN NaN]; %% Leaf Shed Water Stress [0-1]
Tls_L = [NaN NaN NaN NaN NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN NaN NaN NaN NaN NaN]; %%% Stem allocation for Crops  
LAI_min_L = [0.01 0.01 0.01 0.01 0.01 0.01];

Mf_L= [0 0 0 0 0 0]; %% fruit maturation turnover [1/d]
Wm_L= [0 0 0 0 0 0] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.2 0.2 0.2 0.2 0.2 0.2]; %% Allocation to reserve parameter [0-1]
Klf_L =  [1/20 1/20 1/20 1/20 1/20 1/20] ; % [1/83]; %% Dead Leaves fall turnover [1/d]
fab_L = [1.0 1.0 1.0 1.0 1.0 1.0]; %% fraction above-ground sapwood and reserve
fbe_L =[0.0 0.0 0.0 0.0 0.0 0.0]; %% fraction below-ground sapwood and reserve
 
%wheat, barley,grass, Potatoes, RS, Peas.

[ParEx_L(1)]=Exudation_Parameter(0); 
[ParEx_L(2)]=Exudation_Parameter(0); 
[ParEx_L(3)]=Exudation_Parameter(0);
[ParEx_L(4)]=Exudation_Parameter(0);
[ParEx_L(5)]=Exudation_Parameter(0);
[ParEx_L(6)]=Exudation_Parameter(0);
[Mpar_L(1)]=Vegetation_Management_Parameter;
[Mpar_L(2)]=Vegetation_Management_Parameter;
[Mpar_L(3)]=Vegetation_Management_Parameter;
[Mpar_L(4)]=Vegetation_Management_Parameter;
[Mpar_L(5)]=Vegetation_Management_Parameter;
[Mpar_L(6)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%
Mpar_L(1).Date_sowing = [    datenum(2004,01,2,10,0,0)   datenum(2004,9,29,10,0,0) datenum(2005,8,10,10,0,0)  datenum(2006,5,5,10,0,0)   datenum(2006,10,19,10,0,0) datenum(2007,8,28,10,0,0) datenum(2008,10,7,10,0,0) datenum(2009,8,22,10,0,0) datenum(2010,5,1,10,0,0)  datenum(2010,10,15,10,0,0) datenum(2011,9,24,10,0,0) datenum(2012,9,4,10,0,0)  datenum(2013,10,19,10,0,0)  datenum(2014,9,14,10,0,0) datenum(2015,8,5,10,0,0) datenum(2016,5,1,10,0,0) datenum(2016,10,12,10,0,0) datenum(2017,8,30,10,0,0) datenum(2018,10,11,10,0,0)];
Mpar_L(1).Date_harvesting = [datenum(2004,8,4,10,0,0)    datenum(2005,7,14,10,0,0) datenum(2005,11,9,10,0,0)  datenum(2006,10,15,10,0,0) datenum(2007,7,15,10,0,0)  datenum(2008,7,16,10,0,0) datenum(2009,7,21,10,0,0) datenum(2009,11,9,10,0,0) datenum(2010,7,19,10,0,0) datenum(2011,8,2,10,0,0)   datenum(2012,7,9,10,0,0)  datenum(2013,7,28,10,0,0) datenum(2014,7,24,10,0,0)   datenum(2015,7,15,10,0,0) datenum(2015,11,9,10,0,0) datenum(2016,7,25,10,0,0) datenum(2017,7,19,10,0,0) datenum(2018,7,12,10,0,0) datenum(2019,7,19,10,0,0) ];
%%%%%%%%%%%%%%%%%%%%% Winter Wheat
Mpar_L(1).Date_sowing = Mpar_L(1).Date_sowing  ; %% Date of Sowing
Mpar_L(1).Date_harvesting = Mpar_L(1).Date_harvesting; %%% Fully Harvested
Mpar_L(1).Crop_B=[8 10] ;%    %% 1.5 -10 gC m2 
Mpar_L(1).Crop_crown =[1 0 0 0 1 0 1 0 0 1 0 0 1 0 0 0 1 0 1]; 
%%%%%%%%%%  Winter Barley
Mpar_L(2).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(2).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(2).Crop_B=[8 10]; %  
Mpar_L(2).Crop_crown =[0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%%  Grass Cover crops // summer oat, phacelia, Alexandrine clover
Mpar_L(3).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(3).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(3).Crop_B=[8 8]; 
Mpar_L(3).Crop_crown =[0 0 1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%% Potatoes 
Mpar_L(4).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(4).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(4).Crop_B=[10 10]; 
Mpar_L(4).Crop_crown =[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%% Rape Seed
Mpar_L(5).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(5).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(5).Crop_B=[8 12]; % root and then reserves
Mpar_L(5).Crop_crown =[0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0]; 
%%%%%%%%%%%%%
%%%%%%%%%% Peas
Mpar_L(6).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(6).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(6).Crop_B=[8 12];  % second one is C in seed
Mpar_L(6).Crop_crown =[0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0]; 
%%%%%%%%%%%%%
%%%%
Mpar_L(1).fract_left=0;
Mpar_L(1).fract_left_fr=0;
Mpar_L(1).fract_left_AB=0;
Mpar_L(1).fract_left_BG=1;
Mpar_L(2).fract_left=0;
Mpar_L(2).fract_left_fr=0;
Mpar_L(2).fract_left_AB=0;
Mpar_L(2).fract_left_BG=1; 
Mpar_L(3).fract_left=0;
Mpar_L(3).fract_left_fr=0;
Mpar_L(3).fract_left_AB=0;
Mpar_L(3).fract_left_BG=1;
Mpar_L(4).fract_left=0;
Mpar_L(4).fract_left_fr=0;
Mpar_L(4).fract_left_AB=0;
Mpar_L(4).fract_left_BG=1;
Mpar_L(5).fract_left=0;
Mpar_L(5).fract_left_fr=0;
Mpar_L(5).fract_left_AB=0;
Mpar_L(5).fract_left_BG=1;
Mpar_L(6).fract_left=0;
Mpar_L(6).fract_left_fr=0;
Mpar_L(6).fract_left_AB=0;
Mpar_L(6).fract_left_BG=1;

%%%%%%
%%%%%%%%%%%% PRODUCTIVITY
Vmax_H = [0 0 0 0 0 0]; %
%%%%%%%%%%%%%%%%%%%%%%
%%%%
%wheat, barley,grass, Potatoes, RS, Peas.


L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%
%%%%%%%%% DORMANT 1 - MAX GROWTH 2 - NORMAL GROWTH 3 - SENESCENCE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%
LAI_H(1,:)=[0.0 0.0 0.0 0.0 0.0 0.0]; %
B_H(1,:,:)= [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]; %%
Rrootl_H(1,:) = [0 0 0 0 0 0] ;
PHE_S_H(1,:)=[0 0 0 0 0 0];
dflo_H(1,:)=[0 0 0 0 0 0];
AgeL_H(1,:)=[0 0 0 0 0 0];
e_rel_H(1,:)=[0 0 0 0 0 0];
hc_H(1,:) =[0 0 0 0 0 0]; %%
SAI_H(1,:) = [0 0 0 0 0 0]; %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 Live Leaves/ 2 Sapwood/ 3 Fine Roots / 4 Carbohydrate Reserve /5 Fruit and Flower /6 Heartwood - Dead Sapwood / 7 Standing Dead Leaves
LAI_L(1,:)=[0.0 0.0 0.0 0.0 0.0 0.0]; %
B_L(1,:,:)= [0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0]; %%    
Rrootl_L(1,:) = [0 0 0 0 0 0] ;
PHE_S_L(1,:)=[1 1 1 1 1 1]; %[3];
dflo_L(1,:)=[0 0 0 0 0 0];
AgeL_L(1,:)=[0 0 0 0 0 0 ] ;%[730];
e_rel_L(1,:)=[1 1 1 1 1 1];
hc_L(1,:) =[0.15 0.15 0.15 0.15 0.15 0.15]; %%
SAI_L(1,:) = [0.001 0.001 0.001 0.001 0.001 0.001]; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [0 0 0 0 0 0 ];
Preserve_H(1,:)= [0 0 0 0 0 0 ];
Kreserve_H(1,:)= [0 0 0 0 0 0 ];
FNC_H(1,:)=[1 1 1 1 1 1];
NupI_H(1,:,:)= [0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0];
Nreserve_L(1,:)= [1000 1000 1000 1000 1000 1000];
Preserve_L(1,:)= [100 100 100 100 100 100];
Kreserve_L(1,:)= [1000 1000 1000 1000 1000 1000];
FNC_L(1,:)=[1 1 1 1 1 1];
NupI_L(1,:,:)= [0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0 ; 0 0 0];
RexmyI(1,:)= [0 0 0]; 
%%%%%%%%%

if OPT_SoilBiogeochemistry == 1

  
    Nreserve_L(1,:) = 0; %
    Preserve_L(1,:) = 0 ; %
    Kreserve_L(1,:) = 0 ; %
    
    FNC_L(1,:) = 1;
    NupI_L(1,:,:) = repmat([0.100826 0.01101  0.043436947], 6, 1);

    
    RexmyI(1,:,:) = [0.033995 0.1972202  0];
    NavlI(1,:,:) =[  0.244   0.0115  0.0618];% guess the two above are not crop specific ?
    %%%%
    %%%
    %ORIGINAL
%     P(1,:)=   1000*[ 0.124083987036102   0.325573628166407   0.054069975887170                   0                   0   0.101512860630231   0.126564601766440,...
%         0.014062733529605   0.755595303823738   1.096103576915454   6.402023138646436   0.007894696102488   0.013135775003470   0.000162321193831,...
%         0.000104278919687   0.000081160596916   0.000156418379530   0.038738455769914   0.133541727544331   0.065957977030849                   0,...
%         0.002528243856056   0.009999140372626                   0   0.004492687719394   0.913446366042410   0.005555237467077   0.011170245415232,...
%         0.003668130264076                   0   0.000185349798914   0.000017774714467   0.000036321130404   0.000252824385605   0.001098054538178,...
%         0   0.000476491894740   0.137001078180096   0.001107194360735   0.001675533186330   0.000550219539611                   0,...
%         0.000007588306019   0.150000000000000   0.000172964791821   0.015000467306710   0.000000368010905   0.005528056261498                   0,...
%         0.000603437980692   0.016621912220522   0.000015673158772   0.000009163307352   0.000808451498098   0.500387085728242];
%     %%%
    %My test to run from the spin up onwards:
    %could restate here P values from 
%Pini = load("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\TC_4\RUN_RESULTS\RUN_RESULTSCH_OE2_a.mat",'P');
%P(1,:) = Pini.P(15000,:);

Pini = load("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\TCNOV\Outputs\PvariableatendCHOE2April.mat",'P');
P(1,:) = Pini.P(6211,:);


%DepN= 0.00136; %
    %DepP=1.578e-005; %  gP/m2 day
    %DepK=1.27e-004;
    
    %% OLD SHOULD STILL WORK
    %     FertN=0*ones(1,366);
    %     FertN([50 100 200]) = 30;
    %     FertP=0*ones(1,366);
    %     FertP([50 100 200]) = 30;
    %     FertK=0*ones(1,366);
    %     FerKN([50 100 200]) = 30;
    
    %% NEW
    
    Time = datetime( (datenum([1959 1 1]):datenum([2020 1 1]))', "convertfrom", "datenum");
    FertN = timetable(Time,...
        zeros(size(Time)), 'VariableNames', {'Nrate'});
    
    FertP = timetable(Time,...
        zeros(size(Time)), 'VariableNames', {'Prate'});
    
    FertK = timetable(Time,...
        zeros(size(Time)), 'VariableNames', {'Krate'});
    
    %add fertilizer quantities and timing here...

%%%%CAREFUL !!! THE BELOW ONLY WORKS IF WE START POST SPIN UP!!!!!

FertApplication = readtimetable("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\Data\Sites\CH-Oe2(Oensingen crop)\Observations\CH-OE2.Management.allyears.xlsx");
FertDates = datenum(FertApplication.DateYear);
Nrate = FertApplication.N*2;
Prate = FertApplication.P*2;
Krate = FertApplication.K*2; % 

% Ndates = [datenum(2004,3,4)];


%     Ndates = datenum(1959:1:2015,5,20);
%     Nrate  = 50*ones(size(Ndates));
Ndates = [FertDates];
Nrate = [Nrate];
    for nn = 1:length(Ndates)
        ix = find( abs(Ndates(nn) -  datenum(FertN.Time))<0.01);
        if ~isempty(ix)
            FertN.Nrate(ix) = Nrate(nn);
        end
    end    
   

% note - To convert P2O5 to P divide by 2.29.
%     Pdates = datenum(1959:1:2015,5,20);
%     Prate  = 50*ones(size(Pdates));
Pdates = [FertDates];
Prate = [Prate];

    for nn = 1:length(Pdates)
        ix = find( abs(Pdates(nn) -  datenum(FertP.Time))<0.01);
        if ~isempty(ix)
            FertP.Prate(ix) = Prate(nn);
        end
    end
    %note -To convert K2O to K divide by 1.21
%     Kdates = datenum(1959:1:2015,5,20);
%     Krate  = 50*ones(size(Kdates));
    Kdates = [FertDates];
    Krate = [Krate];
    for nn = 1:length(Kdates)
        ix = find( abs(Kdates(nn) -  datenum(FertK.Time))<0.01);
        if ~isempty(ix)
            FertK.Krate(ix) = Krate(nn);
        end

    end

    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1)*1000,Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl);
    B_IO.SC_par=[1 1 1 1];
    
    B_IO.FertN = FertN;
    B_IO.FertP = FertP;
    B_IO.FertK = FertK;
    
    PHs=5;

end
%%%
%%%
TBio_L=[1 1 1 1 1 1];  %%[ton DM / ha ]
TBio_H=[0 0 0 0 0 0];  %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=[0 0 0 0 0 0];  %% [mm/ m2 PFT];
Vl_H=[0 0 0 0 0 0];  %% [mm/ m2 PFT];
Vx_L=[0 0 0 0 0 0];   %% [mm/ m2 PFT];
Vl_L=[1 1 1 1 1 1];   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Initial Conditions
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Ts_under(1)=NaN; 
Tdamp(1)=3;
Tdp(1,:)= 3*ones(1,ms);
TdpI_H(1,:)=3;
TdpI_L(1,:)=3;
%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
In_Litter(1)=0;
SP_wc(1)=0 ; %%[mm]
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
%O(1,:)= [    0.4310    0.4306    0.4304    0.4306    0.4310    0.4317    0.4321    0.4326    0.4327    0.4319    0.4298    0.426,...
    %0.4222    0.4187    0.4138    0.4106    0.4098    0.4105    0.4117    0.4128]; 

O(1,:) = repmat(0.4310, ms, 1);




%dz(ms)=ZWT(1)-Zs(ms);
%%%%%%%%%%%%%%%%%%%%%%%%
%Vmin= 0*dz;
%Vmax= (Osat-Ohy)*dz;
%Vmax(ms) = Vmax(ms) + (Zbed-ZWT(1))*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%V(1,ms)= 0 + (Osat-Ohy)*(Zbed-ZWT(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
cd(cur_dir)
%%%%%%%%%%%%%%%%%
