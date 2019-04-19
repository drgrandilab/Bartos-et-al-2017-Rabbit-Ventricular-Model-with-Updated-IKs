%% Main file.
% This file loads the initial conditions, calls the ode solver, and plots
% the results.
% 
% Reference: D.C. Bartos, S. Morotti, K.S. Ginsburg, E. Grandi, D.M. Bers.
% Quantitative analysis of the Ca2+-dependent regulation of delayed 
% rectifier K+ current IKs in rabbit ventricular myocytes. J Physiol. 2016.
% doi:10.1113/JP273676.
% 
% Please cite the above paper when using this model.

close all
clear all 
clc
%% Parameters for external modules

% ECC and CaM modules
freq = 3%;                  % [Hz] CHANGE DEPENDING ON FREQUENCY
cycleLength = 1e3/freq;     % [ms]
CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM] 
CaNtotDyad = 3e-3/8.293e-4; % [uM] 
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

% Parameters for CaMKII module
% CaMKII expression level: 'WT', 'OE', or 'KO'
expression = 'WT';
CKIIOE = 0; % Should be zero during 'WT' and 'KO' runs
if strcmp(expression,'OE')
    CKIIOE = 1; % Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    CaMKIItotDyad = 120*6;        % [uM] 
    CaMKIItotSL = 120*8.293e-4*6; % [uM]
    CaMKIItotCyt = 120*8.293e-4*6;% [uM]
elseif strcmp(expression,'KO')
    CaMKIItotDyad = 0;          % [uM] 
    CaMKIItotSL = 0;            % [uM]
    CaMKIItotCyt = 0;           % [uM]
end

LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
Ligtot = 0%; % 0.1, 0.05, 0.02 or 0 % [uM] SET ISO CONCENTRATION HERE
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = PLBtot;          % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]
PLMtotBA = 48;              % [uM] - [umol/L cytosol] as in Yang & Saucerman (mouse) model
MyototBA = 70;              % [uM] - [umol/L cytosol] as TnI
IKrtotBA = 0.025;           % [uM] - [umol/L cytosol] as IKs
IClCatotBA = 0.025;         % [uM] - [umol/L cytosol] as ICFTR

% Input parameter for stimulation protocols
prot_input_par = 10; 
%% Collect all parameters and define mass matrix for BAR module

ny_ECC = 83; % 83 state vars in ECC module
ny_cam = 15; % 15 state vars in each CaM module
ny_CaMKII = 6; % 6 state vars in CaMKII module
ny_BAR = 30+6; % 30+6 state vars in BAR module

p = [cycleLength,prot_input_par,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,...
    PLMtotBA,MyototBA,IKrtotBA,IClCatotBA,CKIIOE];

% Need to define Mass matrix for BAR portion
m0 = ones(1,ny_ECC+3*ny_cam+ny_CaMKII); % All state variables in other modules, not BAR
m1 =[0,     0,      0,      1,      1,      1,       1,       1,       1];
m2 =[0,     0,      0,      1,      1,      1];
m3 =[0,     0,      0];
m4 =[1,     1,      0,      0,      1,      1];
m5 =[1,     1,      0,      0,      1,      1,       1,       1      0,      0,      1,      1];

M = diag([m0,m1,m2,m3,m4,m5]); 
%% Establish and define globals

global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store 
global gates_store Jserca_store IKs_store Jleak_store ICFTR_store Incx_store
global I_NaK_store I_kr_store I_Nabk_store
global Lmyo_store Fmyo_store Vmax_store
global IKp_store I_clca_store I_clbk_store Ipca_store Icabk_store Cai_store

tStep = 1; tArray = zeros(1,1e6); I_Ca_store = zeros(1,1e6); I_to_store = zeros(3,1e6);
I_Na_store = zeros(1,1e6); I_K1_store = zeros(1,1e6); ibar_store = zeros(1,1e6);
gates_store = zeros(2,1e6); Jserca_store = zeros(1,1e6); IKs_store = zeros(1,1e6);
Jleak_store = zeros(1e6,2); ICFTR_store = zeros(1,1e6); Incx_store = zeros(1,1e6);
I_NaK_store = zeros(1,1e6); I_kr_store = zeros(1,1e6); I_Nabk_store = zeros(1,1e6);
Fmyo_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6); Vmax_store = zeros(1,1e6);
IKp_store = zeros(1,1e6); I_clca_store = zeros(1,1e6); I_clbk_store = zeros(1,1e6);
Ipca_store = zeros(1,1e6); Icabk_store = zeros(1,1e6); Cai_store = zeros(1,1e6);
%% Load initial conditions

%load yfin_iks_isometric_1Hz % 1-Hz (steady-state), no ISO  
load yfin_iks_isometric_3Hz % 3-Hz (steady-state), no ISO  
%load yfin_iks_isometric_3Hz_0p020iso_240s % 3-Hz, 240 s with 20 nM ISO 
%load yfin_iks_isometric_3Hz_0p050iso_240s % 3-Hz, 240 s with 50 nM ISO 

y0n = yfinal;    
%% Define [Ca] for simulations with fixed [Ca]
% For clamping [Ca] to the initial value, set the flag Ca_clamp = 1 in 
% rabbit_myofilament_masterODEfile and rabbit_myofilament_eccODEfile.

% Ca_conc_clamp = 500*1e-6; % (mM)
% yfinal(38) = Ca_conc_clamp; % yfinal(37) = Ca_conc_clamp; % yfinal(36) = Ca_conc_clamp;
%% Load voltage command for AP-clamp simulations
% For using this signal as voltage command, set protocol = 'AP_clamp' in
% rabbit_myofilament_eccODEfile.

global AP_Em AP_t
load AP_Bartos
%% Run single simulation

tic
tspan = [0 1*1e3]; % [ms]
options = odeset('Mass',M,'RelTol',1e-5,'MaxStep',2); 
[t,y] = ode15s(@rabbit_myofilament_masterODEfile,tspan,y0n,options,p);
yfinal = y(end,:)';
toc
%% Saving final conditions

%save yfin_iks_isometric_1Hz yfinal
%save yfin_iks_isometric_3Hz yfinal
%save yfin_iks_isometric_3Hz_0p020iso_240s yfinal
%save yfin_iks_isometric_3Hz_0p050iso_240s yfinal
%% Rename outputs

tA = tArray(1:tStep); Vmax = Vmax_store(1:tStep);
I_Ca = I_Ca_store(1:tStep); Ito = I_to_store(1,1:tStep);
Itof = I_to_store(2,1:tStep); Itos = I_to_store(3,1:tStep);
INa = I_Na_store(1:tStep); IK1 = I_K1_store(1:tStep);
s1 = gates_store(1,1:tStep); k1 = gates_store(2,1:tStep);
Jserca = Jserca_store(1:tStep); Iks = IKs_store(1:tStep);
Jleak = Jleak_store(1:tStep,:); ICFTR = ICFTR_store(1:tStep);
Incx = Incx_store(1:tStep); INaK = I_NaK_store(1:tStep);
Ikr = I_kr_store(1:tStep); INabk = I_Nabk_store(1:tStep);
Ikp = IKp_store(1:tStep); Iclca = I_clca_store(1:tStep);
Iclbk = I_clbk_store(1:tStep); Ipca = Ipca_store(1:tStep);
Icabk = Icabk_store(1:tStep); Cai = Cai_store(1:tStep);
Lm = Lmyo_store(1:tStep); Fm = Fmyo_store(1:tStep);
%% Plot

color = 'k'; limit_x1 = 0; limit_x2 = 300;
figure, set(gcf,'color','w','Position',[100 100 500 800])
subplot(3,1,1); hold on; set(gca,'box','off','tickdir','out','fontsize',12)
plot(t,y(:,39),color),ylabel('Em (mV)'),xlim([limit_x1 limit_x2])
subplot(3,1,2); hold on; set(gca,'box','off','tickdir','out','fontsize',12)
plot(tA,Iks,color); ylabel('IKs (A/F)'),xlim([limit_x1 limit_x2])
subplot(3,1,3); hold on; set(gca,'box','off','tickdir','out','fontsize',12)
plot(t,y(:,38)*1e6,color),ylabel('[Ca2+]i (nM)'),xlim([limit_x1 limit_x2])
xlabel('Time (ms)');