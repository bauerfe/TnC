%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MAIN PARAMATERS WE CHANGE  
%[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]
Sl_L = [0.028 0.025 0.023 0.04 0.024 0.027]; % SLA [m^2 gC] specific leaf area of  biomass [m^2 /gC]
r_L = [0.025 0.023 0.027 0.025 0.020 0.026];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.19 0.28 0.31 0.25 0.24 0.28]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
Tcold_L = [3 0 3 -2 3 3]; %% [°C] Cold Leaf SLed
dmg_L= [65 75 110 60 50 50 ]; %%% Tree 30 Grasses Day of Max Growth %[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]
mjDay_L = [200 222 200 350 250 250]; %% Maximum Julian day for leaf onset
LDay_min_L =[12.5 12.3 12.3 10.2 11.1 10.2]; %% Minimum Day duration for leaf onset WW ok 9.7,
age_cr_L= [95 90 45 85 120 90];%[120 125 250 90 120 90]; %% [day] Critical Leaf Age %[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]
Tlo_L = [13.5 2 13.4 9 17 3.5]; %% Mean Temperature for Leaf onset WW ok with 8 
LDay_cr_L = [5 7.7 7 7 9 10]; %%%  Threshold for senescence day light [h] %%%%%% 
Sl_emecrop_L= [0.01 0.022 0.021 0.011 0.015 0.001]; %%% Additional SLA at emergence  [ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]
MHcrop_L =[1.2 1.2 0.7 1.2 3 1.2]; %%[m] maximum crop height  
Vmax_L = [88 106 65 35 45 54]; % %[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]
LtR_L = [1.0 1.1 1.2 1.0 1.1 1.6]; %%% Leaf to Root ratio maximum
Trr_L = [8.5 7.5 6.9 6.5 8.5 6.5]; %% Translocation rate [gC /m^2 d]

%In order to play with yields, the parameter soCrop_L, and ff_r_L are good targets

soCrop_L = [0.10 0.10 0.0 0.05 0.05 0.05]; %
ff_r_L= [0.50 0.4 0.50 0.50 0.52 0.50]; % fraction allocated to fruit - to the grain pool 5


% SLA used in CLM SB and pot was 0.2, WW was 0.28

% Sugar Beet	30/03/2004	29/09/2004 (1)
% Winter Wheat	14/10/2004	03/08/2005 (2)
% Potatoes	01/05/2006	19/09/2006 (3)
% Winter Wheat	13/10/2006	05/08/2007 (2)
% Sugar Beet	22/04/2008	04/11/2008 (1)
% Winter Wheat	13/11/2008	07/08/2009 (2)
% Mustard	01/09/2009	01/12/2009 (4)
% Potatoes	24/04/2010	03/09/2010 (3)
% Winter Wheat	13/10/2010	16/08/2011 (2)
% Maize	14/05/2012	13/10/2012 (5)
% Winter Wheat	25/10/2012	12/08/2013 (2)
% Mustard	05/09/2013	15/11/2013 (4)
% Potatoes	05/04/2014	21/08/2014 (3)
% Winter Wheat	14/10/2014	02/08/2015 (2)
% Mustard	21/08/2015	07/12/2015 (4)
% Sugar Beet	12/04/2016	27/10/2016 (1)
% Winter Wheat	29/10/2016	29/07/2017 (2)
% Mustard	06/09/2017	06/12/2017 (4)
% Potatoes	22/04/2018	11/09/2018 (3) % defoliant which we dont use...
% Winter Wheat	10/10/2018	01/08/2019 (2)
% Oat/Faba Bean	09/08/2019	04/12/2019 (6)
% Sugar Beet	01/04/2020	12/11/2020 (1)

% Sugar Beet (1)
% Winter Wheat (2)
% Potatoes (3)
% Mustard (4)
% Maize (5)
% Oat (6)

%[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]

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
Cbare =  1; Ccrown = [0.0 0.0 0.0 0.0 0.0 0.0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL INPUT %%%%  silty clay loam
Pcla= 0.18; % actually 18.2 
Psan= 0.05; % actually 5.10 and then 68% silt. - check this... 
Porg= 0.10; %0.05 between 5 and 10 % 
Color_Class = 0; % silt is 70%, sand is 5%, clay is 20%

% Pcla= 0.40;
% Psan= 0.10;
% Porg= 0.025;
% Color_Class = 0; % copied from CHOE2


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
Zs = [0 10 20 30 80 90 1250 1500 2500 4000 5000]; %%% [ms+1] %% need to ensure this aligns with OBS SWC readings depths
Zs = [0 100 200 500 800 1000 1250 1500 2500 4000 5000]; %%% [ms+1] % other results which more or less work...

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
ZR95_L = [950 950 950 950 950 950]; %% [mm] 
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
FI_L=[0.081 0.081 0.081 0.081 0.081 0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons] % note should be 0.04 for C4 which is maize...
Do_L=[1000 1000 1000 1000 1000 1000]; %%[Pa]
a1_L=[8 8 8 8 6 8]; % Empirical parameter connecting stomatal aperture and net assimilation
go_L=[0.01 0.01 0.01 0.01 0.005 0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3 3 3 3 4 3 ];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants C3 or C4
DSE_L =[0.649 0.649 0.649 0.649 0.649 0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72 72 72 72 72 72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf Inf Inf Inf Inf Inf];
rjv_L= [2.4 2.4 2.4 2.4 2.4 2.4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pss_H = [800]; Pwp_H = [3000]; %%% [kPa]
%Pss_L = [500]; Pwp_L = [3500]; %%% [kPa]
Psi_sto_50_H =  [-2.0 -2.0 -2.0 -2.0 -2.5 -2.0] ;%% [MPa]  Water Potential at 50% loss conductivity
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
Psi_sto_00_L = [-0.7 -0.7 -0.7 -0.7 -0.5 -0.7];%  %% [MPa]  Water Potential at PLCs loss conductivity % increasing doesnt fix it decreasing doesnt either
Psi_sto_50_L = [-3 -3  -2  -3  -3 -3 ]-1;%  %% [MPa]  Water Potential at 50% loss conductivity [-2.9 -2.9  -2.9  -2.9  -2.9  -2.9 ]
%%% Leaf
PsiL00_L =  [-1.2 -1.2 -1.2 -1.2 -1.2 -1.2]; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L =  [-5.5 -5.5 -5.5 -5.5 -5.5 -5.5];%%[MPa]  Water Potential at 50% loss conductivity  [-4.5 -5.5 -4.5 -4.5 -4.5 -4.5]
Kleaf_max_L = [5 5 5 5 5 5 ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200 1200 1200 1200 1200];  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [0.0 0.0 0.0 0.0 0.0 0.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000 80000 80000 80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-5.5 -5.5 -5.5 -5.5 -0.01 -5.5]; %%[MPa]  Water Potential at 50% loss conductivity [-4.5 -4.5 -4.5 -4.5 -4.5 -4.5]
Cx_L= [150 150 150 150 150 150]; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 

%%%% Growth Parameters
PsiG50_H= [-2.45 -2.45 -2.45 -2.45 -2.45 -2.45];  %%[MPa]
PsiG99_H= [-4.2 -4.2 -4.2 -4.2 -4.2 -4.2];  %%[MPa]
gcoef_H = [3.5 3.5 3.5 3.5 3.5 3.5]; % [gC/m2 day]
%%------

PsiG50_L = [-1.2 -1.2 -1.2 -1.2 -1.2 -1.2]-1;
PsiG99_L= [-3.5 -3.5 -3.5 -3.5 -3.5 -3.5]-1;
gcoef_L = [3.5 3.5 3.5 3.5 3.5 3.5]; % [gC/m2 day]

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
Sl_H = [0.016 0.016 0.016 0.016 0.016 0.016]; % 0.018 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
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



%%


Nl_L= [16 16 16 16 16 16]; %[gC/gN ] Leaf Carbon-Nitrogen ratio % had a look where to find this, tricky, so left at 16
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
dc_C_L = [7/365 3/365 7/365 7/365 7/365 7/365]; %% [1/ d°C] -- [Factor of increasing mortality]
drn_L=  [1/365 1/365 1/365 1/365 1/365 1/365]; %% turnover root  [1/d]
dsn_L= [0 0 0 0 0 0]; % normal transfer rate sapwood [1/d]


Bfac_lo_L= [0.99 0.99 0.99 0.99 0.99 0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN NaN NaN NaN NaN NaN]; %% Leaf Shed Water Stress [0-1]


Tls_L = [NaN NaN NaN NaN NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN NaN NaN NaN NaN NaN]; %%% Stem allocation for Crops  


%[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]


LAI_min_L = [0.01 0.01 0.01 0.01 0.01 0.01];

Mf_L= [0 0 0 0 0 0]; %% fruit maturation turnover [1/d]
Wm_L= [0 0 0 0 0 0] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.2 0.2 0.2 0.2 0.2 0.2]; %% Allocation to reserve parameter [0-1]
%[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]

Klf_L =  [1/20 1/20 1/20 1/20 1/20 1/20] ; % [1/83]; %% Dead Leaves fall turnover [1/d]
fab_L = [1.0 1.0 1.0 1.0 1.0 1.0]; %% fraction above-ground sapwood and reserve
fbe_L =[0.0 0.0 0.0 0.0 0.0 0.0]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.50 0.50 0.50 0.50 0.52 0.50]; % fraction allocated to fruit - to the grain pool 5 I presume it is 

%[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat]




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

% % this has been inputed here on 18-02-24 - could be the issue...

% Mpar_L(3).LAI_cut =[0.04]; %% LAI of grass after cut %% 2.65 Height of 15cm more or less - was 1.68
% Mpar_L(3).jDay_cut = [datenum([2006 7 4])];
% 
% Mpar_L(3).LAI_cut =[0.04]; %% LAI of grass after cut %% 2.65 Height of 15cm more or less - was 1.68
% Mpar_L(3).jDay_cut = [datenum([2010 7 4])];
% 
% Mpar_L(3).LAI_cut =[0.04]; %% LAI of grass after cut %% 2.65 Height of 15cm more or less - was 1.68
% Mpar_L(3).jDay_cut = [datenum([2014 7 4])];
% 
% Mpar_L(3).LAI_cut =[0.04]; %% LAI of grass after cut %% 2.65 Height of 15cm more or less - was 1.68
% Mpar_L(3).jDay_cut = [datenum([2018 7 4])];


% Mpar_L(1).Date_sowing = [ datenum(2004,01,2,10,0,0)  datenum(2004,9,29,10,0,0) datenum(2005,8,19,10,0,0) datenum(2006,5,5,10,0,0) datenum(2006,10,19,10,0,0) datenum(2007,8,28,10,0,0) datenum(2008,10,7,10,0,0) datenum(2009,8,22,10,0,0) datenum(2010,5,9,10,0,0) datenum(2010,10,15,10,0,0) datenum(2011,9,24,10,0,0) datenum(2012,9,4,10,0,0) datenum(2013,10,19,10,0,0)  datenum(2014,9,14,10,0,0) datenum(2015,8,15,10,0,0) datenum(2016,5,9,10,0,0) datenum(2016,10,12,10,0,0) datenum(2017,8,30,10,0,0) datenum(2018,10,11,10,0,0)];
% Mpar_L(1).Date_harvesting = [datenum(2004,8,4,10,0,0) datenum(2005,7,14,10,0,0) datenum(2005,11,9,10,0,0)  datenum(2006,10,15,10,0,0) datenum(2007,7,15,10,0,0) datenum(2008,7,16,10,0,0) datenum(2009,7,21,10,0,0) datenum(2009,11,9,10,0,0) datenum(2010,7,19,10,0,0) datenum(2011,8,2,10,0,0) datenum(2012,7,9,10,0,0) datenum(2013,7,28,10,0,0) datenum(2014,7,24,10,0,0) datenum(2015,7,15,10,0,0) datenum(2015,11,9,10,0,0) datenum(2016,7,25,10,0,0) datenum(2017,7,19,10,0,0) datenum(2018,7,12,10,0,0) datenum(2019,7,19,10,0,0) ];

% Mpar_L(1).Date_sowing = [ datenum(2004,4,8,10,0,0) datenum(2004,10,14,10,0,0) datenum(2006,5,15,10,0,0) datenum(2006,10,13,10,0,0) datenum(2008,4,29,10,0,0) datenum(2008,11,13,10,0,0) datenum(2009,9,15,10,0,0) datenum(2010,5,8,10,0,0) datenum(2010,10,13,10,0,0) datenum(2012,6,20,10,0,0) datenum(2012,10,25,10,0,0) datenum(2013,9,15,10,0,0) datenum(2014,5,9,10,0,0) datenum(2014,10,14,10,0,0) datenum(2015,08,21,10,0,0) datenum(2016,04,18,10,0,0) datenum(2016,10,29,10,0,0) datenum(2017,09,6,10,0,0) datenum(2018,04,29,10,0,0) datenum(2018,10,10,10,0,0) datenum(2019,08,9,10,0,0) datenum(2020,04,1,10,0,0) ];
% Mpar_L(1).Date_harvesting = [datenum(2004,9,29,10,0,0) datenum(2005,8,3,10,0,0) datenum(2006,9,19,10,0,0) datenum(2007,8,5,10,0,0) datenum(2008,11,4,10,0,0)  datenum(2009,8,7,10,0,0) datenum(2009,12,1,10,0,0) datenum(2010,9,3,10,0,0)  datenum(2011,8,11,10,0,0) datenum(2012,10,13,10,0,0) datenum(2013,8,12,10,0,0) datenum(2013,11,15,10,0,0) datenum(2014,8,14,10,0,0) datenum(2015,8,2,10,0,0) datenum(2015,12,7,10,0,0) datenum(2016,10,27,10,0,0) datenum(2017,07,29,10,0,0) datenum(2017,12,6,10,0,0) datenum(2018,09,11,10,0,0) datenum(2019,08,01,10,0,0) datenum(2019,12,04,10,0,0) datenum(2020,11,12,10,0,0) ];
%  
%real sowing and harvest are below
Mpar_L(1).Date_sowing =     [datenum(2004,3,30,10,0,0) datenum(2004,10,14,10,0,0)   datenum(2006,5,1,10,0,0)  datenum(2006,10,13,10,0,0)   datenum(2008,4,22,10,0,0)   datenum(2008,11,13,10,0,0)   datenum(2009,9,1,10,0,0) datenum(2010,4,29,10,0,0) datenum(2010,10,13,10,0,0) datenum(2012,5,14,10,0,0) datenum(2012,10,25,10,0,0) datenum(2013,9,15,10,0,0) datenum(2014,5,9,10,0,0) datenum(2014,10,14,10,0,0) datenum(2015,08,21,10,0,0) datenum(2016,04,12,10,0,0) datenum(2016,10,29,10,0,0) datenum(2017,09,6,10,0,0) datenum(2018,04,22,10,0,0) datenum(2018,10,10,10,0,0) datenum(2019,08,9,10,0,0) datenum(2020,04,1,10,0,0) ];
Mpar_L(1).Date_harvesting = [datenum(2004,9,29,10,0,0) datenum(2005,08,03,10,0,0)   datenum(2006,9,19,10,0,0) datenum(2007,08,05,10,0,0)   datenum(2008,11,4,10,0,0)   datenum(2009,08,07,10,0,0)   datenum(2009,12,1,10,0,0) datenum(2010,9,3,10,0,0)  datenum(2011,8,11,10,0,0) datenum(2012,10,13,10,0,0) datenum(2013,8,12,10,0,0) datenum(2013,11,15,10,0,0) datenum(2014,8,14,10,0,0) datenum(2015,8,2,10,0,0) datenum(2015,12,7,10,0,0) datenum(2016,10,27,10,0,0) datenum(2017,07,29,10,0,0) datenum(2017,12,6,10,0,0) datenum(2018,09,11,10,0,0) datenum(2019,08,01,10,0,0) datenum(2019,12,04,10,0,0) datenum(2020,11,12,10,0,0) ];
 
Mpar_L(1).Date_harvesting = [datenum(2004,9,29,10,0,0) datenum(2005,08,03,10,0,0)   datenum(2006,8,3,10,0,0) datenum(2007,08,05,10,0,0)   datenum(2008,11,4,10,0,0)   datenum(2009,08,07,10,0,0)   datenum(2009,12,1,10,0,0) datenum(2010,8,15,10,0,0)  datenum(2011,8,11,10,0,0) datenum(2012,10,13,10,0,0) datenum(2013,8,12,10,0,0) datenum(2013,11,15,10,0,0) datenum(2014,7,5,10,0,0) datenum(2015,8,2,10,0,0) datenum(2015,12,7,10,0,0) datenum(2016,10,27,10,0,0) datenum(2017,07,29,10,0,0) datenum(2017,12,6,10,0,0) datenum(2018,07,15,10,0,0) datenum(2019,08,01,10,0,0) datenum(2019,12,04,10,0,0) datenum(2020,11,12,10,0,0) ];
% above is with potato defoliant added, trying to see how to best replicate
% this...


Mpar_L(1).Date_sowing =     [datenum(2004,3,30,10,0,0) datenum(2004,10,14,10,0,0)   datenum(2006,5,1,10,0,0)  datenum(2006,10,13,10,0,0)   datenum(2008,4,22,10,0,0)   datenum(2008,11,13,10,0,0)   datenum(2009,9,1,10,0,0) datenum(2010,4,29,10,0,0) datenum(2010,10,20,10,0,0) datenum(2012,5,14,10,0,0) datenum(2012,10,25,10,0,0) datenum(2013,9,22,10,0,0) datenum(2014,5,9,10,0,0) datenum(2014,10,14,10,0,0) datenum(2015,08,29,10,0,0) datenum(2016,04,12,10,0,0) datenum(2016,10,29,10,0,0) datenum(2017,09,13,10,0,0) datenum(2018,04,22,10,0,0) datenum(2018,10,10,10,0,0) datenum(2019,08,9,10,0,0) datenum(2020,04,1,10,0,0) ];




% 2006 is aug 3
%Aug 15 2010
%Julz 21 2014
%July 22 2018 % for defoliant for potatoes



% Sugar Beet	30/03/2004	29/09/2004 (1) 1 
% Winter Wheat	14/10/2004	03/08/2005 (2) 2
% Potatoes	01/05/2006	19/09/2006 (3) 3
% Winter Wheat	13/10/2006	05/08/2007 (2) 4
% Sugar Beet	22/04/2008	04/11/2008 (1) 5
% Winter Wheat	13/11/2008	07/08/2009 (2)6
% Mustard	01/09/2009	01/12/2009 (4) 7
% Potatoes	24/04/2010	03/09/2010 (3) 8
% Winter Wheat	13/10/2010	16/08/2011 (2) 9
% Maize	14/05/2012	13/10/2012 (5) 10
% Winter Wheat	25/10/2012	12/08/2013 (2) 11
% Mustard	05/09/2013	15/11/2013 (4) 12
% Potatoes	05/04/2014	21/08/2014 (3) 13
% Winter Wheat	14/10/2014	02/08/2015 (2) 14
% Mustard	21/08/2015	07/12/2015 (4) 15 15
% Sugar Beet	12/04/2016	27/10/2016 (1) 16
% Winter Wheat	29/10/2016	29/07/2017 (2) 17
% Mustard	06/09/2017	06/12/2017 (4) 18 - this is where it blocks
% Potatoes	22/04/2018	11/09/2018 (3) 19
% Winter Wheat	10/10/2018	01/08/2019 (2) 20
% Oat/Faba Bean	09/08/2019	04/12/2019 (6) 21
% Sugar Beet	01/04/2020	12/11/2020 (1) 22

%Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat
%%%%%%%%%%%%%%%%%%%%% Sugar Beet 1 
Mpar_L(1).Date_sowing = Mpar_L(1).Date_sowing  ; %% Date of Sowing
Mpar_L(1).Date_harvesting = Mpar_L(1).Date_harvesting; %%% Fully Harvested
Mpar_L(1).Crop_B=[8 10] ;%    %% 1.5 -10 gC m2 
Mpar_L(1).Crop_crown =[1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1]; %22 crops
%%%%%%%%%%  Winter Wheat 2
Mpar_L(2).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(2).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(2).Crop_B=[8 20]; %  
Mpar_L(2).Crop_crown =[0 1 0 1 0 1 0 0 1 0 1 0 0 1 0 0 1 0 0 1 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%%  Potatoes 3
Mpar_L(3).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(3).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(3).Crop_B=[8 8]; 
Mpar_L(3).Crop_crown =[0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%% Mustard 4
Mpar_L(4).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(4).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(4).Crop_B=[10 10]; 
Mpar_L(4).Crop_crown =[0 0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%% Maize 5 
Mpar_L(5).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(5).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(5).Crop_B=[8 10]; 
Mpar_L(5).Crop_crown =[0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]; 
%%%%%%%%%%%%%
%%%%%%%%%% Oat
Mpar_L(6).Date_sowing = Mpar_L(1).Date_sowing ; %% Date of Sowing
Mpar_L(6).Date_harvesting = Mpar_L(1).Date_harvesting ; %%% Fully Harvested
Mpar_L(6).Crop_B=[8 10]; 
Mpar_L(6).Crop_crown =[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]; 
%%%%%%%%%%%%%
%%%%
Mpar_L(1).fract_left=0;
Mpar_L(1).fract_left_fr=0;
Mpar_L(1).fract_left_AB=0;
Mpar_L(1).fract_left_BG=1; % only leave the roots... 

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


%Vmax_L = [95 110 90 100 45 61]; % %[ Sugar Beet, Winter Wheat, Potatoes, Mustard, Maize, Oat] - wheat was 125 which is unrealistic 

%%%%%%%%%%%%%%%%%%%%%%
%%%%

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
%     %ORIGINAL
%     P(1,:)=   1000*[ 0.124083987036102   0.325573628166407   0.054069975887170                   0                   0   0.101512860630231   0.126564601766440,...
%         0.014062733529605   0.755595303823738   1.096103576915454   6.402023138646436   0.007894696102488   0.013135775003470   0.000162321193831,...
%         0.000104278919687   0.000081160596916   0.000156418379530   0.038738455769914   0.133541727544331   0.065957977030849                   0,...
%         0.002528243856056   0.009999140372626                   0   0.004492687719394   0.913446366042410   0.005555237467077   0.011170245415232,...
%         0.003668130264076                   0   0.000185349798914   0.000017774714467   0.000036321130404   0.000252824385605   0.001098054538178,...
%         0   0.000476491894740   0.137001078180096   0.001107194360735   0.001675533186330   0.000550219539611                   0,...
%         0.000007588306019   0.150000000000000   0.000172964791821   0.015000467306710   0.000000368010905   0.005528056261498                   0,...
%         0.000603437980692   0.016621912220522   0.000015673158772   0.000009163307352   0.000808451498098   0.500387085728242];
%     %%%
%     My test to run from the spin up onwards:
%     could restate here P values from 

% Pini = load('Pvariable.mat', 'P');
% P(1,:) = Pini.P(6209,:);


% Pini = load("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\TC_4\RUN_RESULTS\RUN_RESULTSCH_OE2_a.mat",'P');
% P(1,:) = Pini.P(15000,:);

% Pini = load("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\TCNOV\Outputs\Pvariableatend.mat");
% P(1,:) = Pini.P(6210,:);

Pini = load('PvariableSpinUPBELON.mat', 'P');
P(1,:) = Pini.P(6210,:);

% 
% Pini = load('Pvariable.mat', 'P');
% P(1,:) = Pini.P(1,:);


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

%      Before the period under study,
% 2003: 180 kg ha1 mineral N and 60 kg ha1 organic N (sugar lime)
% 2002: 156 kg ha1 mineral N, 14 kg ha1
% phosphorus and 42 kg ha1 potassium in 2002, 
% 2001: 180 kg ha1 mineral N in 2001. (Aubinet et al. 2009)

%%%%CAREFUL !!! THE BELOW ONLY WORKS IF WE START POST SPIN UP!!!!!

FertApplication = readtimetable("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\Data\Sites\BE-Lon(Lonzee)\Inputs\FertBELON.xlsx");

FertDates = datenum(FertApplication.DateYear);
Nrate = FertApplication.N;
Prate = FertApplication.P;
Krate = FertApplication.K;

%"In total, 201 kg mineral nitrogen per hectare was applied in four
%fractions."
% Ndates = [datenum(2004,3,4)];


%     Ndates = datenum(1959:1:2015,5,20);
%     Nrate  = 50*ones(size(Ndates));
Ndates = [FertDates];
Nrate = [Nrate];
    for nn = 1:length(Ndates) % 43
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
    %%
    
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
