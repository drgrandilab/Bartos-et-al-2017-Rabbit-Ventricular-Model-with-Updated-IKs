function ydot = rabbit_myofilament_eccODEfile(t,y,p)
% This function describes the excitation-contraction coupling.

%% State variables
% 1       2       3       4       5       6       7       8      9       10  
% m       h       j       xks_sl  f       fcaBj   fcaBsl  xtos   ytos    xtof  
% 11      12      13      14      15      16      17      18     19      20   
% ytof    xkr     xks_j   RyRr    RyRo    RyRi    NaBj    NaBsl  TnCL    TnCHc
% 21      22      23      24      25      26      27      28     29      30
% TnCHm   CaM     Myoc    Myom    SRB     SLLj    SLLsl   SLHj   SLHsl   Csqnb
% 31      32      33      34      35      36      37      38     39      40
% Ca_sr   Naj     Nasl    Nai     Ki      Caj     Casl    Cai    Vm      rtos
% 41: INaL h gate
% 42-53: placeholder for INa Markov model
% 54-59: myofilament (TSCa; TSCa_star; TSCa_tilde; TS_star; L_p; L_w)
% 60-83: ICaL Markov model - see lines 530-537

ydot = zeros(size(y));
%% Assign passed-in parameters

cycleLength = p(1);
% CaMKII phosphorylated targets (%)
LCC_CKp = p(2);
RyR_CKp = p(3);
PLB_CKp = p(4);
% PKA phosphorylated targets (%)
LCCa_PKAp = p(5);
LCCb_PKAp = p(6);
PLB_PKAn = p(7); % non-phosphorylated PLB targets
RyR_PKAp = p(8);
TnI_PKAp = p(9);
IKs_PKAp = p(10);
ICFTR_PKAp = p(11);
PLM_PKAp = p(12);
Myo_PKAp = p(13);
IKr_PKAp = p(14);
IClCa_PKAp = p(15);
% Flag for CKII OE
CKIIOE = p(16);
% Protocol parameter
input_par = p(17);
%% Define Simulation Properties

% Protocol selection
protocol = 'pace'; % current-clamp
%protocol = 'v_step'; % voltage-clamp
%protocol = 'AP_clamp'; % AP-clamp

% Set Ca_clamp = 1 for clamping [Ca] to the initial value, 0 otherwise.
% (define accordingly Ca_clamp in masterODEfile)
Ca_clamp = 0;
%% Other flags

% Set CKIIflag to 1 for CKII OE, 0 otherwise
CKIIflag = CKIIOE;
% Set ICa_MarkovFlag to 1 for Markov ICa, 0 otherwise
ICa_MarkovFlag = 1;
% Set MarkovFlag to 1 for Markov INa, 0 otherwise
MarkovFlag = 0;
% Set Ito to use either original params (=0) or Grandi Params (=1)
ItoFlag = 1;
% Set mechFlag to 0 for isometric or 1 for isotonic contraction
mechFlag = 0;
% Set IKs_flag to 0 for original IKs model
IKs_flag = 1;
%% Model Parameters

% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K] 310 K (37 C) for BT / 295 K (22 C) for RT
FoRT = Frdy/R/Temp;
Cmem = 1.3810e-10;   % [F] membrane capacitance
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um]
cellRadius = 10.25;   % cell radius [um]
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]
distSLcyto = 0.45;    % dist. SL to cytosol [um]
distJuncSL = 0.5;  % dist. junc to SL [um]
DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec] 
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]
J_ca_juncsl = 1/1.2134e12; % [L/msec]
J_ca_slmyo = 1/2.68510e11; % [L/msec]
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] 
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] 

% Fractional currents in compartments
Fjunc = 0.11; Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]
%% Na transport parameters

GNa =  16.0;        % [mS/uF]
GNaB = 0.297e-3;    % [mS/uF] 
IbarNaK = 1.90719;     % [uA/uF]
KmNaip = 11;         % [mM]
KmKo = 1.5;         % [mM]
Q10NaK = 1.63;  
Q10KmNai = 1.39;
%% K currents parameters

pNaK = 0.01833;      
GtoSlow = 0.06;     % [mS/uF] 
GtoFast = 0.02;     % [mS/uF] 
gkp = 0.001;
%% Cl current parameters

GClCa = 0.109625;   % [mS/uF]
GClB = 9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]
%% Ca transport parameters

pNa = 1.5e-8;       % [cm/sec]
pCa = 5.4e-4;       % [cm/sec]
pK = 2.7e-7;        % [cm/sec]
KmCa = 0.6e-3;      % [mM]
Q10CaL = 1.8;       

IbarNCX = 9.0;      % [uA/uF]
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact = 0.256e-3;   % [mM] 
Q10NCX = 1.57;      % [none]
IbarSLCaP = 0.0673; % [uA/uF]
KmPCa = 0.5e-3;     % [mM] 
GCaB = 2.513e-4;    % [uA/uF] 
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 2.86e-4;  % [mM/msec] (mmol/L cytosol/msec)
Kmf = 0.246e-3;          % [mM]
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10;               % [mM^-2 1/ms]      
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]
%% Buffering parameters

Bmax_Naj = 7.561;       % [mM] 
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = 19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] 
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.38e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.62e-3*Vmyo/Vjunc*0.1;    % [mM]    
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.35e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] 
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 2.7;        %140e-3*Vmyo/Vsr; [mM] 
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 
%% Myofilament parameters

nc=3;                 %Addim.(=1,2,3...)
Ap=1008e+04;          %[mN/mm2/um/mM]
Aw=Ap/5;              %[mN/mm2/um/mM]
alfa=0.5;             %[mN/mm2]
bet=80;               %[1/um]
Bp=0.5;               %[1/ms]
Bw=0.35;              %[1/ms]
f=0.0023;             %[1/ms]
gama=28000;           %[1/um2]
hpr=0.006;            %[um]
hwr=0.0001;           %[um]
Ke=105000;            %[mN/mm2/um^5] 
Lz=0.97;              %[um] (0.97)
La=1.15;              %[um]
Lc=1.05;              %[um]
Le=10;                %[mN/mm2/um](10)
RLa=20;               %[1/um2]
TSt=0.07/nc;          %[mM]
Yb=181.6e+06;
Yc=4;                 %[1/um]
Yd=0.028;             %[1/ms]
Yp=0.1397;            %[1/ms]  
Yr=0.1397;            %[1/ms]
Yv=0.9;               %[1/ms]
Za=0.0023;            %[1/ms]
Zb=0.1397;            %[1/ms]
Zp=0.2095;            %[1/ms]
Zr=7262.6e+06;        %[1/mM^nc/ms]
Fh=0.1;               %Addim.
%% Global Variable for Time

global tStep tArray               
if t > tArray(tStep)    % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;  
end
tArray(tStep) = t;  
%% I_Na: Fast Na Current

% Max INa alterations with CKII hyperactivity as in Hund & Rudy 2008
if CKIIflag == 1    
    inashift = -3.25;
    alphaCKII = -.18;
    deltGbarNal_CKII = 2;  
else
    inashift = 0;
    alphaCKII = 0;
    deltGbarNal_CKII = 0;
end

am = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bm = 0.08*exp(-y(39)/11);
if (y(39)-inashift) >= -40
    ah = 0; aj = 0;
    bh = 1/(0.13*(1+exp(-((y(39) - inashift)+10.66)/11.1)));
    bj = 0.3*exp(-2.535e-7*(y(39)-inashift))/(1+exp(-0.1*((y(39)-inashift)+32)));
else
    ah = 0.135*exp((80+(y(39)-inashift))/-6.8);
    bh = 3.56*exp(0.079*(y(39)-inashift))+3.1e5*exp(0.35*(y(39)-inashift));
    % Including alteration to aj as in Hund and Rudy 2008
    aj = (1+alphaCKII)*((-1.2714e5*exp(0.2444*(y(39)-inashift))-3.474e-5*exp(-0.04391*(y(39)-inashift)))*((y(39)-inashift)+37.78)/(1+exp(0.311*((y(39)-inashift)+79.23))));
    bj = 0.1212*exp(-0.01052*(y(39)-inashift))/(1+exp(-0.1378*((y(39)-inashift)+40.14)));
end
ydot(1) = am*(1-y(1))-bm*y(1);
ydot(2) = ah*(1-y(2))-bh*y(2);
ydot(3) = aj*(1-y(3))-bj*y(3);
I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc);
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);
I_Na = I_Na_junc+I_Na_sl;
%% I_Na,L: Late Na current (as in Hund & Rudy 2008)

GbarNal = .0065*(1+deltGbarNal_CKII)*2;   % deltGbar assigned in 'Fast INa' section

% h-gate (note - m gate is same as INa m gate - using y(1) for this)
hlss = 1/(1+exp((y(39)+91)/6.1));
tauhl = 600; % ms
ydot(41) = (hlss-y(41))/tauhl;
I_Nalj = Fjunc*GbarNal*y(1)^3*y(41)*(y(39)-ena_junc);
I_Nalsl = Fsl*GbarNal*y(1)^3*y(41)*(y(39)-ena_sl);
I_Nal = I_Nalj+I_Nalsl;

% Reassign total junctional and sl currents as sum of fast and late currents
I_Na_junc = I_Na_junc+I_Nalj;
I_Na_sl = I_Na_sl+I_Nalsl;
%% Placeholders for Markov INa 

% If MarkovFlag = 1, INa currents will be dictated by Markov scheme. If
% flag = 0, the original H-H scheme is used to compute the current.

% Parameters
GNa2 = 23;      % [mS/uF]
% If not using INa Markov, set odes to zero to speed simulations
% ydot(42:53) = zeros(1,12); % already set to 0, see line 45

I_Na_junc2 = Fjunc*GNa2*(y(55)+y(59))*(y(39)-ena_junc);     % junctional current
I_Na_sl2 = Fsl*GNa2*(y(55)+y(59))*(y(39)-ena_sl);           % sl current
%% Compute total INa (fast and late components)

I_Na_junc = (1-MarkovFlag)*I_Na_junc + MarkovFlag*I_Na_junc2;
I_Na_sl = (1-MarkovFlag)*I_Na_sl + MarkovFlag*I_Na_sl2;

I_Na = I_Na_junc + I_Na_sl;

global I_Na_store
I_Na_store(tStep) = I_Na;
%% I_nabk: Na Background Current

I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

global I_Nabk_store
I_Nabk_store(tStep) = I_nabk;
%% I_nak: Na/K Pump Current
%PLM_PKAp = 0.006294;

% PKA-dependent PLM phosphoregulation
fracPKA_PLMo = 0.006294; % Derived quantity (PLM_PKAp(baseline)/PLMtot)
fracPKA_PLMiso = 0.9232; % Derived quantity (PLM_PKAp(ISO)/PLMtot)
kPKA_PLM=KmNaip*(1-13.6/18.8)/(fracPKA_PLMiso/fracPKA_PLMo-1); % PLM_PKAp ISO
    %kPKA_PLM=KmNaip*(0.5)*(1-13.6/18.8)/(fracPKA_PLMiso/fracPKA_PLMo-1); % PLM_PKAp ISO (50%)
KmNaip_PKA=-kPKA_PLM+kPKA_PLM*(PLM_PKAp/fracPKA_PLMo);
KmNaip = KmNaip-KmNaip_PKA; % 27.66% reduction w/ ISO

sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

global I_NaK_store
I_NaK_store(tStep) = I_nak;
%% I_kr: Rapidly Activating K Current

% PKA-dependent IKr phosphoregulation
fracPKA_Ikro = 0.1098; % Derived quantity (IKr_PKAp(baseline)/Ikrtot)
fracPKA_Ikriso = 0.8380; % Derived quantity (IKr_PKAp(ISO)/Ikrtot)
kPKA_Ikr=(IKr_PKAp-fracPKA_Ikro)/(fracPKA_Ikriso-fracPKA_Ikro);
Vkr_PKA = 10*kPKA_Ikr;% 10 mV shift w/ ISO
dGkr_PKA = 0.3*kPKA_Ikr;% 30% increase w/ ISO 
% total effect: SSA shifted and peak increased by 37%
%Vkr_PKA = 0; dGkr_PKA = 0;

gkr = (1+dGkr_PKA)*0.03*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+50+Vkr_PKA)/7.5));
tauxr = 1/(1.38e-3*(y(39)+7)/(1-exp(-0.123*(y(39)+7)))+6.1e-4*(y(39)+10)/(exp(0.145*(y(39)+10))-1));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+33)/22.4));
I_kr = gkr*y(12)*rkr*(y(39)-ek);

global I_kr_store
I_kr_store(tStep) = I_kr;
%% I_ks: Slowly Activating K Current

% PKA-dependent IKs phosphoregulation % 100 nM ISO: IKs_PKAp=0.8380;
fracPKA_Ikso = 0.1098; % Derived quantity (IKs_PKAp(baseline)/Ikrtot)
fracPKA_Iksiso = 0.8380; % Derived quantity (IKs_PKAp(ISO)/Ikrtot)
kPKA_Iks=(IKs_PKAp-fracPKA_Ikso)/(fracPKA_Iksiso-fracPKA_Ikso);

if IKs_flag == 0,% OLD IKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PKA-dependent IKs phosphoregulation (as specified in Negroni et al, table 2)
    fracIKsavail = 1+(2.8-1)*kPKA_Iks; % 2.8-fold increase w/ ISO
    Xs05 = 1.5-(13.5+1.5)*kPKA_Iks; % 15 mV shift w/ ISO

    pcaks_junc = -log10(y(36))+3.0; 
    pcaks_sl = -log10(y(37))+3.0;  
    gks_junc = fracIKsavail*0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6))); % Now regulated by PKA
    gks_sl = fracIKsavail*0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));     % Now regulated by PKA
    eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));	
    % xsss = 1/(1+exp(-(y(39)-1.5)/16.7)); % Original version
    xsss = 1/(1+exp(-(y(39)-Xs05)/16.7)); % Now regulated by PKA
    tauxs = 1/(7.19e-5*(y(39)+30)/(1-exp(-0.148*(y(39)+30)))+1.31e-4*(y(39)+30)/(exp(0.0687*(y(39)+30))-1)); 
    ydot(13) = (xsss-y(13))/tauxs;
    I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
    I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);
else % NEW IKS - 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PKA-dependent IKs phosphoregulation
    gks_factor = 0.05;
    P_g_0 = gks_factor*(0.2+0.2*kPKA_Iks); % 0.2->0.4 w/ ISO
    P_g_max = gks_factor*(0.8+0.5*kPKA_Iks); % 0.8->1.3 w/ ISO
    P_vh_0 = -1-10*kPKA_Iks; % -1->-11 w/ ISO
    P_vh_max = -12-9*kPKA_Iks; % -12->-21 w/ ISO
    P_tau_0 = 26+9*kPKA_Iks; % 160->260 w/ ISO
    P_tau_max = 40+4*kPKA_Iks; % 300->370 w/ ISO

    caks_junc = y(36); caks_sl = y(37); % normal simulation
        
    gks_junc = P_g_0 + (P_g_max-P_g_0)/(1 + (150e-6/caks_junc)^1.3); % Regulated by PKA
    gks_sl = P_g_0 + (P_g_max-P_g_0)/(1 + (150e-6/caks_sl)^1.3); % Regulated by PKA
    VsXs_Ca_junc = P_vh_0 + (P_vh_max-P_vh_0)/(1 + (350e-6/caks_junc)^4); % Regulated by PKA
    xsss_junc = 1/(1+exp(-(y(39)-VsXs_Ca_junc)/25));
    VsTs_Ca_junc = P_tau_0 + (P_tau_max-P_tau_0)/(1 + (150e-6/caks_junc)^3); % Regulated by PKA
    tauxs_junc = 2*(50+(50+350*exp(-((y(39)+30)^2)/4000))*1/(1+exp(-(y(39)+VsTs_Ca_junc)/10)));
    VsXs_Ca_sl = P_vh_0 + (P_vh_max-P_vh_0)/(1 + (350e-6/caks_sl)^4); % Regulated by PKA
    xsss_sl = 1/(1+exp(-(y(39)-VsXs_Ca_sl)/25));
    VsTs_Ca_sl = P_tau_0 + (P_tau_max-P_tau_0)/(1 + (150e-6/caks_sl)^3); % Regulated by PKA
    tauxs_sl = 2*(50+(50+350*exp(-((y(39)+30)^2)/4000))*1/(1+exp(-(y(39)+VsTs_Ca_sl)/10))); 
    ydot(13) = (xsss_junc-y(13))/tauxs_junc;
    ydot(4) = (xsss_sl-y(4))/tauxs_sl;
    eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));
    I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
    I_ks_sl = Fsl*gks_sl*y(4)^2*(y(39)-eks);
end

I_ks = I_ks_junc+I_ks_sl;

global IKs_store
IKs_store(tStep) = I_ks;
%% I_kp: Plateau K current

kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

global IKp_store
IKp_store(tStep) = I_kp;
%% I_to: Transient Outward K Current (slow and fast components)

xtoss = 1/(1+exp(-(y(39)+3.0)/15));
ytoss = 1/(1+exp((y(39)+33.5)/10));
rtoss = 1/(1+exp((y(39)+33.5)/10));
tauxtos = 9/(1+exp((y(39)+3.0)/15))+0.5;

if ItoFlag == 0
% Shannon Versions
    tauytos = 3e3/(1+exp((y(39)+60.0)/10))+30;
    taurtos = 2.8e3/(1+exp((y(39)+60.0)/10))+220; 
elseif ItoFlag == 1 && CKIIflag == 0
    % Grandi Versions
    Py = 182; Pr1 = 8085; Pr2 = 313;                % Normal
    tauytos = Py/(1+exp((y(39)+33.5)/10))+1;
    taurtos = Pr1/(1+exp((y(39)+33.5)/10))+Pr2;     % Should Pr2 be 313 or 131 (as in Grandi paper)?
elseif ItoFlag == 1 && CKIIflag == 1
    Py = 15; Pr1 = 3600; Pr2 = 500; GtoSlow = GtoSlow*1.5;  % CKII OE
    tauytos = Py/(1+exp((y(39)+33.5)/10))+1;
    taurtos = Pr1/(1+exp((y(39)+33.5)/10))+Pr2;   
end

ydot(8) = (xtoss-y(8))/tauxtos;
ydot(9) = (ytoss-y(9))/tauytos;
ydot(40)= (rtoss-y(40))/taurtos; 
I_tos = GtoSlow*y(8)*(y(9)+0.5*y(40))*(y(39)-ek);    % [uA/uF]
% tauxtof = 3.5*exp(-y(39)*y(39)/30/30)+1.5;    % Original
tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5;       % Version in Grandi Code (does not change AP shape)
tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = GtoFast*y(10)*y(11)*(y(39)-ek);
I_to = I_tos + I_tof;

global I_to_store
I_to_store(1,tStep) = I_to;     % Total I_to
I_to_store(2,tStep) = I_tof;    % Fast Component
I_to_store(3,tStep) = I_tos;    % Slow component
%% I_k1: Time-Independent K Current

aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);
I_ki = 0.9*sqrt(Ko/5.4)*kiss*(y(39)-ek);

global I_K1_store
I_K1_store(tStep) = I_ki;
%% I_ClCa & I_Clbk: Ca-activated Cl Current & background Cl Current

% PKA-dependent IClCa phosphoregulation
fracPKA_IClCao = 0.1624; % Derived quantity (IClCa_PKAp(baseline)/Ikrtot)
fracPKA_IClCaiso = 0.8918; % Derived quantity (IClCa_PKAp(ISO)/Ikrtot)
kPKA_IClCa=(IClCa_PKAp-fracPKA_IClCao)/(fracPKA_IClCaiso-fracPKA_IClCao);
KdClCa = (1+(0.704-1)*kPKA_IClCa)*KdClCa; % -29.6% w/ 100 nM ISO

I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;

I_Clbk = GClB*(y(39)-ecl);

global I_clca_store I_clbk_store
I_clca_store(tStep) = I_ClCa;     % Total I_ClCa
I_clbk_store(tStep) = I_Clbk;     % I_Clbk
%% I_Ca: Original H-H formulation for LTCC (unused if ICa_MarkovFlag = 1)

% If not using H&H model, set odes to zero to speed simulations
dss = 1/(1+exp(-(y(39)+14.5)/6.0));
taud = dss*(1-exp(-(y(39)+14.5)/6.0))/(0.035*(y(39)+14.5));
fss = 1/(1+exp((y(39)+35.06)/3.6))+0.6/(1+exp((50-y(39))/20));
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+14.5))^2 )+0.02);
% ydot(4) = (dss-y(4))/taud;
% ydot(5) = (fss-y(5))/(tauf);
% ydot(6) = (1.7)*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc  
% ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl
% ydot(4:7) = zeros(1,4); % already set to 0, see line 45
% NOW y(4) is used for IKs !!!
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);

I_Ca_junc1 = (Fjunc_CaL*ibarca_j*y(4)*y(5)*(1-y(6))*Q10CaL^Qpow)*0.45;
I_Ca_sl1 = (Fsl_CaL*ibarca_sl*y(4)*y(5)*(1-y(7))*Q10CaL^Qpow)*0.45;
I_CaK1 = (ibark*y(4)*y(5)*(Fjunc_CaL*(1-y(6))+Fsl_CaL*(1-y(7)))*Q10CaL^Qpow)*0.45;
I_CaNa_junc1 = (Fjunc_CaL*ibarna_j*y(4)*y(5)*(1-y(6))*Q10CaL^Qpow)*0.45;
I_CaNa_sl1 = (Fsl_CaL*ibarna_sl*y(4)*y(5)*(1-y(7))*Q10CaL^Qpow)*.45;
%% I_Ca: LTCC MARKOV MODEL - based on Mahajan et al. (2008)

% This portion contains Markov state transitions for four channel types:
% 'mode 1' channels in the junction and sl and 'mode 2' channels in the
% same two compartments. Markov state transitions are computed for each
% channel type independently - total currents are the sum of the two
% channel types in each compartment (i.e. ICatot_junc = ICa_mode1_junc +
% ICa_mode2_junc). Ca-dependent transition rates differ between juncitonal
% and sl channels, whereas closing rate (r2) is adjusted to define mode1
% vs. mode2 channels. Parameters determined through microscopic
% reversibility are redefined to preserve constraint.

% CaMKII shifts distribution of junctional and subsarcolemmal channels to 
% either mode 2 at the expense of mode 1 channels (i.e. 10% mode 2 results 
% in 90% mode 1).

% PKA alters overall availability of channels (favail term that changes
% overall scaling factor for currents) and also shifts distribution of
% mode1/2 channels. PKA actions act on both junctional and sarcolemmal
% channels.

% input LTCC module ODE (24 state vars)
Pc2_LCCj_m1=y(60); Pc1_LCCj_m1=y(61); Pi1Ca_LCCj_m1=y(62);
Pi2Ca_LCCj_m1=y(63); Pi1Ba_LCCj_m1=y(64); Pi2Ba_LCCj_m1=y(65);
Pc2_LCCj_m2=y(66); Pc1_LCCj_m2=y(67); Pi1Ca_LCCj_m2=y(68);
Pi2Ca_LCCj_m2=y(69); Pi1Ba_LCCj_m2=y(70); Pi2Ba_LCCj_m2=y(71);
Pc2_LCCsl_m1=y(72); Pc1_LCCsl_m1=y(73); Pi1Ca_LCCsl_m1=y(74);
Pi2Ca_LCCsl_m1=y(75); Pi1Ba_LCCsl_m1=y(76); Pi2Ba_LCCsl_m1=y(77);
Pc2_LCCsl_m2=y(78); Pc1_LCCsl_m2=y(79); Pi1Ca_LCCsl_m2=y(80);
Pi2Ca_LCCsl_m2=y(81); Pi1Ba_LCCsl_m2=y(82); Pi2Ba_LCCsl_m2=y(83);

Po_LCCj_m1 =  1-Pc2_LCCj_m1-Pc1_LCCj_m1-Pi1Ca_LCCj_m1-Pi2Ca_LCCj_m1-Pi1Ba_LCCj_m1-Pi2Ba_LCCj_m1;
Po_LCCj_m2 =  1-Pc2_LCCj_m2-Pc1_LCCj_m2-Pi1Ca_LCCj_m2-Pi2Ca_LCCj_m2-Pi1Ba_LCCj_m2-Pi2Ba_LCCj_m2;
Po_LCCsl_m1 = 1-Pc2_LCCsl_m1-Pc1_LCCsl_m1-Pi1Ca_LCCsl_m1-Pi2Ca_LCCsl_m1-Pi1Ba_LCCsl_m1-Pi2Ba_LCCsl_m1;
Po_LCCsl_m2 = 1-Pc2_LCCsl_m2-Pc1_LCCsl_m2-Pi1Ca_LCCsl_m2-Pi2Ca_LCCsl_m2-Pi1Ba_LCCsl_m2-Pi2Ba_LCCsl_m2;

% To allow for CDI KO
cajLCC = y(36);
caslLCC = y(37);

% LCC Current Fixed Parameters
pCa = 5.4e-4;       % [cm/sec] - Ca permeability
taupo = 1;          % [ms] - Time constant of activation
TBa = 450;          % [ms] - Time constant
s1o = .0221;
k1o = .03;
kop = 2.5e-3;       % [mM]
cpbar = 8e-3;       % [mM]
tca = 78.0312;
ICa_scale = 5.25;
recoveryReduc = 3;

%%% PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph) %%%%%%
%fracLCCbpo = .0328; % Derived quantity - (LCCbp(baseline)/LCCbtot)
%favail = 1*(.017*LCCb_PKAp/fracLCCbpo + 0.983);   % Test (max x1.5 pCa)
fracLCCbpo = .03284; % Derived quantity - (LCCbp(baseline)/LCCbtot)
favail = 1+0.5*(LCCb_PKAp - fracLCCbpo)/(1 - fracLCCbpo); % 1.5 with max phosphorylation
ICa_scale =  ICa_scale*favail;

% Voltage- and Ca-dependent Parameters
poss = 1/(1+exp(-y(39)/8));
fcaj = 1/(1+(kop/cajLCC)^3);            
Rv = 10 + 4954*exp(y(39)/15.6);
PrLCC = 1-1/(1+exp(-(y(39)+40)/4));     
PsLCC = 1/(1+exp(-(y(39)+40)/11.32));
TCaj = (tca + 0.1*(1+(cajLCC/cpbar)^2))/(1+(cajLCC/cpbar)^2); 
tauCaj = (Rv-TCaj)*PrLCC + TCaj;     
tauBa = (Rv-TBa)*PrLCC + TBa;

% Tranisition Rates (20 rates)
alphaLCC = poss/taupo;
betaLCC = (1-poss)/taupo;
r1 = 0.3;                               % [1/ms] - Opening rate
r2 = 3;                                 % [1/ms] - closing rate
s1 = s1o*fcaj; 
s1p = .00195;                           % [ms] - Inactivation rate
k1 = k1o*fcaj;  
k1p = .00413;                           % [ms] - Inactivation rate
k2 = 1e-4;                              % [ms] - Inactivation rate
k2p = .00224;                           % [ms] - Inactivation rate
s2 = s1*(k2/k1)*(r1/r2);
s2p = s1p*(k2p/k1p)*(r1/r2);
k3 = exp(-(y(39)+40)/3)/(3*(1+exp(-(y(39)+40)/3)));
k3p = k3;
k5 = (1-PsLCC)/tauCaj;
k6 = (fcaj*PsLCC)/tauCaj;
k5p = (1-PsLCC)/tauBa;

% Recovery terms
k5 = k5/recoveryReduc;
k5p = k5p/recoveryReduc;

k6p = PsLCC/tauBa;
k4 = k3*(alphaLCC/betaLCC)*(k1/k2)*(k5/k6);
k4p = k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p);

global gates_store
gates_store(1,tStep) = s1;
gates_store(2,tStep) = k1;

% State transitions for MODE 1 junctional LCCs %%%
% O = no differential; C2 = 60; C1 = 61; I1Ca = 62; I2Ca = 63; I1Ba = 64; I2Ba = 65;
dPc2_LCCj_m1 = betaLCC*Pc1_LCCj_m1 + k5*Pi2Ca_LCCj_m1 + k5p*Pi2Ba_LCCj_m1 - (k6+k6p+alphaLCC)*Pc2_LCCj_m1;                      % C2_m1j
dPc1_LCCj_m1 = alphaLCC*Pc2_LCCj_m1 + k2*Pi1Ca_LCCj_m1 + k2p*Pi1Ba_LCCj_m1 + r2*Po_LCCj_m1 - (r1+betaLCC+k1+k1p)*Pc1_LCCj_m1;   % C1_m1j
dPi1Ca_LCCj_m1 = k1*Pc1_LCCj_m1 + k4*Pi2Ca_LCCj_m1 + s1*Po_LCCj_m1 - (k2+k3+s2)*Pi1Ca_LCCj_m1;                              % I1Ca_m1j
dPi2Ca_LCCj_m1 = k3*Pi1Ca_LCCj_m1 + k6*Pc2_LCCj_m1 - (k4+k5)*Pi2Ca_LCCj_m1;                                                 % I2Ca_m1j
dPi1Ba_LCCj_m1 = k1p*Pc1_LCCj_m1 + k4p*Pi2Ba_LCCj_m1 + s1p*Po_LCCj_m1 - (k2p+k3p+s2p)*Pi1Ba_LCCj_m1;                        % I1Ba_m1j
dPi2Ba_LCCj_m1 = k3p*Pi1Ba_LCCj_m1 + k6p*Pc2_LCCj_m1 - (k5p+k4p)*Pi2Ba_LCCj_m1;                                             % I2Ba_m1j

ibarca_jm1 = (4*pCa*y(39)*Frdy*FoRT)*(.001*exp(2*y(39)*FoRT)-0.341*Cao)/(exp(2*y(39)*FoRT)-1);
I_Ca_junc_m1 = (Fjunc_CaL*ibarca_jm1*Po_LCCj_m1*Q10CaL^Qpow)*ICa_scale;

%%% Re-define all parameters as mode 2 specific parameters %%%
s1om2 = .0221;
k1om2 = .03;
kopm2 = 2.5e-3;
cpbarm2 = 8e-3;
tcam2 = 78.0312;

possm2 = 1/(1+exp(-y(39)/8));
fcajm2 = 1/(1+(kopm2/cajLCC)^3);    % Depends on junctional Ca
Rvm2 = 10 + 4954*exp(y(39)/15.6);
PrLCCm2 = 1-1/(1+exp(-(y(39)+40)/4));     % Correct version I believe
PsLCCm2 = 1/(1+exp(-(y(39)+40)/11.32));
TCajm2 = (tcam2 + 0.1*(1+(cajLCC/cpbarm2)^2))/(1+(cajLCC/cpbarm2)^2); % Caj dependent
tauCajm2 = (Rvm2-TCajm2)*PrLCCm2 + TCajm2;     % Caj dependence
tauBam2 = (Rvm2-TBa)*PrLCCm2 + TBa;

alphaLCCm2 = possm2/taupo;
betaLCCm2 = (1-possm2)/taupo;
r1m2 = 0.3;                               % [1/ms] - Opening rate
r2m2 = 3/10;                                 % [1/ms] - closing rate
s1m2 = s1om2*fcajm2; 
s1pm2 = .00195;                           % [ms] - Inactivation rate
k1m2 = k1om2*fcajm2; 
k1pm2 = .00413;                           % [ms] - Inactivation rate
k2m2 = 1e-4;                              % [ms] - Inactivation rate
k2pm2 = .00224;                           % [ms] - Inactivation rate
s2m2 = s1m2*(k2m2/k1m2)*(r1m2/r2m2);
s2pm2 = s1pm2*(k2pm2/k1pm2)*(r1m2/r2m2);
k3m2 = exp(-(y(39)+40)/3)/(3*(1+exp(-(y(39)+40)/3)));
k3pm2 = k3m2;
k5m2 = (1-PsLCCm2)/tauCajm2;
k6m2 = (fcajm2*PsLCCm2)/tauCajm2;
k5pm2 = (1-PsLCCm2)/tauBam2;
k5m2 = k5m2/recoveryReduc;      % reduced for recovery
k5pm2 = k5pm2/recoveryReduc;    % reduced for recovery    
k6pm2 = PsLCCm2/tauBam2;
k4m2 = k3m2*(alphaLCCm2/betaLCCm2)*(k1m2/k2m2)*(k5m2/k6m2);
k4pm2 = k3pm2*(alphaLCCm2/betaLCCm2)*(k1pm2/k2pm2)*(k5pm2/k6pm2);

%%% State transitions for MODE 2 junctional LCCs %%%
% O = no differential; C2 = 66; C1 = 67; I1Ca = 68; I2Ca = 69; I1Ba = 70; I2Ba = 71;
dPc2_LCCj_m2 = betaLCCm2*Pc1_LCCj_m2 + k5m2*Pi2Ca_LCCj_m2 + k5pm2*Pi2Ba_LCCj_m2 - (k6m2+k6pm2+alphaLCCm2)*Pc2_LCCj_m2;                          % C2_m2j
dPc1_LCCj_m2 = alphaLCCm2*Pc2_LCCj_m2 + k2m2*Pi1Ca_LCCj_m2 + k2pm2*Pi1Ba_LCCj_m2 + r2m2*Po_LCCj_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*Pc1_LCCj_m2;   % C1_m2j
dPi1Ca_LCCj_m2 = k1m2*Pc1_LCCj_m2 + k4m2*Pi2Ca_LCCj_m2 + s1m2*Po_LCCj_m2 - (k2m2+k3m2+s2m2)*Pi1Ca_LCCj_m2;                                  % I1Ca_m2j
dPi2Ca_LCCj_m2 = k3m2*Pi1Ca_LCCj_m2 + k6m2*Pc2_LCCj_m2 - (k4m2+k5m2)*Pi2Ca_LCCj_m2;                                                         % I2Ca_m2j
dPi1Ba_LCCj_m2 = k1pm2*Pc1_LCCj_m2 + k4pm2*Pi2Ba_LCCj_m2 + s1pm2*Po_LCCj_m2 - (k2pm2+k3pm2+s2pm2)*Pi1Ba_LCCj_m2;                            % I1Ba_m2j
dPi2Ba_LCCj_m2 = k3pm2*Pi1Ba_LCCj_m2 + k6pm2*Pc2_LCCj_m2 - (k5pm2+k4pm2)*Pi2Ba_LCCj_m2;                                                     % I2Ba_m2j

ibarca_jm2 = (4*pCa*y(39)*Frdy*FoRT)*(.001*exp(2*y(39)*FoRT)-0.341*Cao)/(exp(2*y(39)*FoRT)-1);
I_Ca_junc_m2 = (Fjunc_CaL*ibarca_jm2*(Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale;

%%% CaMKII AND PKA-DEPENDENT SHIFTING OF DYADIC LCCS TO MODE 2 %%%%
%fpkam2 = 0.1543*LCCa_PKAp - .0043;  % Assumes max phosphorylation results in 15% mode 2 channels
fracLCCapo = .02784; % Derived quantity - (LCCap(baseline)/LCCatot)
fpkam2 = 0.15*(LCCa_PKAp - fracLCCapo)/(1 - fracLCCapo); % Assumes max phosphorylation results in 15% mode 2 channels
if fpkam2 < 0,
    fpkam2 = 0;
end
fckiim2 = LCC_CKp*.1; % Assumes max phosphorylation results in 10% mode 2 channels
% Sum up total fraction of CKII and PKA-shifted mode 2 channels
junc_mode2 = fckiim2 + fpkam2;
% Total junctional ICa
I_Ca_junc2 = (1-junc_mode2)*I_Ca_junc_m1 + junc_mode2*I_Ca_junc_m2;

%%% SUB-SARCOLEMMAL LCCs %%%

% Re-assign necessary params to be Casl sensitive
fcasl = 1/(1+(kop/caslLCC)^3);    % Depends on sl Ca
TCasl = (tca + 0.1*(1+(caslLCC/cpbar))^2)/(1+(caslLCC/cpbar)^2);
tauCasl = (Rv-TCasl)*PrLCC + TCasl;

% Re-assign necessary rates to be Casl sensitive
s1sl = s1o*fcasl;
k1sl = k1o*fcasl;
s2sl = s1sl*(k2/k1sl)*(r1/r2);
s2psl = s1p*(k2p/k1p)*(r1/r2);
k5sl = (1-PsLCC)/tauCasl;
k5sl = k5sl/recoveryReduc;  % Reduced for recovery
k6sl = (fcasl*PsLCC)/tauCasl;
k4sl = k3*(alphaLCC/betaLCC)*(k1sl/k2)*(k5sl/k6sl);
k4psl = k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p);

% State transitions for 'mode 1' sarcolemmal LCCs
% O = no differential; C2 = 72; C1 = 73; I1Ca = 74; I2Ca = 75; I1Ba = 76; I2Ba = 77;
dPc2_LCCsl_m1 = betaLCC*Pc1_LCCsl_m1 + k5sl*Pi2Ca_LCCsl_m1 + k5p*Pi2Ba_LCCsl_m1 - (k6sl+k6p+alphaLCC)*Pc2_LCCsl_m1;                      % C2_m1sl
dPc1_LCCsl_m1 = alphaLCC*Pc2_LCCsl_m1 + k2*Pi1Ca_LCCsl_m1 + k2p*Pi1Ba_LCCsl_m1 + r2*Po_LCCsl_m1 - (r1+betaLCC+k1sl+k1p)*Pc1_LCCsl_m1;    % C1_m1sl
dPi1Ca_LCCsl_m1 = k1sl*Pc1_LCCsl_m1 + k4sl*Pi2Ca_LCCsl_m1 + s1sl*Po_LCCsl_m1 - (k2+k3+s2sl)*Pi1Ca_LCCsl_m1;                         % I1Ca_m1sl
dPi2Ca_LCCsl_m1 = k3*Pi1Ca_LCCsl_m1 + k6sl*Pc2_LCCsl_m1 - (k4sl+k5sl)*Pi2Ca_LCCsl_m1;                                               % I2Ca_m1sl
dPi1Ba_LCCsl_m1 = k1p*Pc1_LCCsl_m1 + k4psl*Pi2Ba_LCCsl_m1 + s1p*Po_LCCsl_m1 - (k2p+k3p+s2psl)*Pi1Ba_LCCsl_m1;                       % I1Ba_m1sl
dPi2Ba_LCCsl_m1 = k3p*Pi1Ba_LCCsl_m1 + k6p*Pc2_LCCsl_m1 - (k5p+k4psl)*Pi2Ba_LCCsl_m1;                                               % I2Ba_m1sl

ibarca_slm1 = (4*pCa*y(39)*Frdy*FoRT)*(.001*exp(2*y(39)*FoRT)-0.341*Cao)/(exp(2*y(39)*FoRT)-1);
I_Casl_m1 = (Fsl_CaL*ibarca_slm1*Po_LCCsl_m1*Q10CaL^Qpow)*ICa_scale;

% Adjust closing rate for 'mode 2' sarcolemmal LCCs
r2slm2 = r2m2;
s2slm2 = s1sl*(k2/k1sl)*(r1/r2slm2);
s2pslm2 = s1p*(k2p/k1p)*(r1/r2slm2);

%%% State transitions for mode 2 sarcolemmal LCCs
% O = no differential; C2 = 78; C1 = 79; I1Ca = 80; I2Ca = 81; I1Ba = 82; I2Ba = 83
dPc2_LCCsl_m2 = betaLCC*Pc1_LCCsl_m2 + k5sl*Pi2Ca_LCCsl_m2 + k5p*Pi2Ba_LCCsl_m2 - (k6sl+k6p+alphaLCC)*Pc2_LCCsl_m2;                      % C2_m2sl
dPc1_LCCsl_m2 = alphaLCC*Pc2_LCCsl_m2 + k2*Pi1Ca_LCCsl_m2 + k2p*Pi1Ba_LCCsl_m2 + r2slm2*Po_LCCsl_m2 - (r1+betaLCC+k1sl+k1p)*Pc1_LCCsl_m2;% C1_m2sl
dPi1Ca_LCCsl_m2 = k1sl*Pc1_LCCsl_m2 + k4sl*Pi2Ca_LCCsl_m2 + s1sl*Po_LCCsl_m2 - (k2+k3+s2slm2)*Pi1Ca_LCCsl_m2;                       % I1Ca_m2sl
dPi2Ca_LCCsl_m2 = k3*Pi1Ca_LCCsl_m2 + k6sl*Pc2_LCCsl_m2 - (k4sl+k5sl)*Pi2Ca_LCCsl_m2;                                               % I2Ca_m2sl
dPi1Ba_LCCsl_m2 = k1p*Pc1_LCCsl_m2 + k4psl*Pi2Ba_LCCsl_m2 + s1p*Po_LCCsl_m2 - (k2p+k3p+s2pslm2)*Pi1Ba_LCCsl_m2;                     % I1Ba_m2sl
dPi2Ba_LCCsl_m2 = k3p*Pi1Ba_LCCsl_m2 + k6p*Pc2_LCCsl_m2 - (k5p+k4psl)*Pi2Ba_LCCsl_m2;                                               % I2Ba_m2sl

ibarca_slm2 = (4*pCa*y(39)*Frdy*FoRT)*(.001*exp(2*y(39)*FoRT)-0.341*Cao)/(exp(2*y(39)*FoRT)-1);
I_Casl_m2 = (Fsl_CaL*ibarca_slm2*Po_LCCsl_m2*Q10CaL^Qpow)*ICa_scale;

% Sum mode 1 and mode 2 sl channels for total sl current
fckiim2_sl = 0; % Set to zero since SL LCCp by CaMKII is negligible
sl_mode2 = fckiim2_sl + fpkam2;
I_Ca_sl2 = (1-sl_mode2)*I_Casl_m1 + sl_mode2*I_Casl_m2; 

% Na and K currents through LCC
I_CaKj2 = ibark*Fjunc_CaL*((1-junc_mode2)*Po_LCCj_m1 + junc_mode2*Po_LCCj_m2)*Q10CaL^Qpow*ICa_scale; 
I_CaKsl2 = ibark*Fsl_CaL*((1-sl_mode2)*Po_LCCsl_m1 + sl_mode2*Po_LCCsl_m2)*Q10CaL^Qpow*ICa_scale;
I_CaK2 = I_CaKj2+I_CaKsl2;
I_CaNa_junc2 = (Fjunc_CaL*ibarna_j*((1-junc_mode2)*Po_LCCj_m1+junc_mode2*Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale;
I_CaNa_sl2 = Fsl_CaL*ibarna_sl*((1-sl_mode2)*Po_LCCsl_m1 + sl_mode2*Po_LCCsl_m2)*Q10CaL^Qpow*ICa_scale;

% These are now able to switch depending on whether or not the flag to
% switch to Markov model of ICa is ON
I_Ca_junc = (1-ICa_MarkovFlag)*I_Ca_junc1 + ICa_MarkovFlag*I_Ca_junc2;
I_Ca_sl = (1-ICa_MarkovFlag)*I_Ca_sl1 + ICa_MarkovFlag*I_Ca_sl2;
I_Ca = I_Ca_junc+I_Ca_sl;   % Total Ca curren throuhgh LCC
I_CaNa_junc = (1-ICa_MarkovFlag)*(I_CaNa_junc1) + (ICa_MarkovFlag)*(I_CaNa_junc2);
I_CaNa_sl = (1-ICa_MarkovFlag)*(I_CaNa_sl1) + (ICa_MarkovFlag)*(I_CaNa_sl2);
I_CaNa = I_CaNa_junc + I_CaNa_sl;   % Total Na current through LCC
I_CaK = (1-ICa_MarkovFlag)*(I_CaK1) + ICa_MarkovFlag*(I_CaK2);  % Total K current through LCC

% Collect all currents through LCC
I_Catot = I_Ca+I_CaK+I_CaNa;

% output LTCC module ODE
ydot(60:65)=[dPc2_LCCj_m1 dPc1_LCCj_m1 dPi1Ca_LCCj_m1 dPi2Ca_LCCj_m1 dPi1Ba_LCCj_m1 dPi2Ba_LCCj_m1];
ydot(66:71)=[dPc2_LCCj_m2 dPc1_LCCj_m2 dPi1Ca_LCCj_m2 dPi2Ca_LCCj_m2 dPi1Ba_LCCj_m2 dPi2Ba_LCCj_m2];
ydot(72:77)=[dPc2_LCCsl_m1 dPc1_LCCsl_m1 dPi1Ca_LCCsl_m1 dPi2Ca_LCCsl_m1 dPi1Ba_LCCsl_m1 dPi2Ba_LCCsl_m1];
ydot(78:83)=[dPc2_LCCsl_m2 dPc1_LCCsl_m2 dPi1Ca_LCCsl_m2 dPi2Ca_LCCsl_m2 dPi1Ba_LCCsl_m2 dPi2Ba_LCCsl_m2];

global I_Ca_store ibar_store
I_Ca_store(tStep) = I_Catot;
ibar_store(tStep) = ibarca_j;
%% I_ncx: Na/Ca Exchanger flux

Ka_junc = 1/(1+(Kdact/y(36))^3);
Ka_sl = 1/(1+(Kdact/y(37))^3);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = (KmCai*Nao^3*(1+(y(32)/KmNai)^3)+KmNao^3*y(36)+ KmNai^3*Cao*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36))*(1+ksat*exp((nu-1)*y(39)*FoRT));
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = (KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)+KmNai^3*Cao*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37))*(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc;
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl;
I_ncx = I_ncx_junc+I_ncx_sl;

global Incx_store
Incx_store(tStep) = I_ncx;
%% I_pca: Sarcolemmal Ca Pump Current

I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

global Ipca_store
Ipca_store(tStep) = I_pca;
%% I_cabk: Ca Background Current

I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

global Icabk_store
Icabk_store(tStep) = I_cabk;
%% I_CFTR or I_cl_(cAMP): Cystic Fibrosis Transmembrane Conductance Reg.

% This is an Em- and time-independent current that is activated by PKA
%fact_pka_cftr = 1.1933*ICFTR_PKAp - 0.1933;
fracCFTRpo = .1624; % Derived quantity - (CFTRp(baseline)/CFTRtot)
fact_pka_cftr = 1*(ICFTR_PKAp - fracCFTRpo)/(1 - fracCFTRpo);
if fact_pka_cftr < 0,
    fact_pka_cftr = 0;
end
gCFTR = fact_pka_cftr*4.9e-3; % [A/F] - Max value as in Shannon et al. (2005) with max phosphorylation
Icftr = gCFTR*(y(39) - ecl);

global ICFTR_store
ICFTR_store(tStep) = Icftr;
%% RyR model - SR Ca release fluxes and leak

MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
kleak = 5.348e-6;

%%% CaMKII and PKA-dependent phosphoregulation of RyR Po %%%
fCKII_RyR = (20*RyR_CKp/3 - 1/3);
%fPKA_RyR = RyR_PKAp*1.025 + 0.9750;
fracRyRpo = .02443; % Derived quantity - (RyRp(baseline)/RyRtot)
fPKA_RyR = 1 + (RyR_PKAp-fracRyRpo) / (1-fracRyRpo); % 2 with max phosphorylation
koSRCa = (fCKII_RyR + fPKA_RyR - 1)*koSRCa;

% ODEs for RyR states and SR release through open RyRs
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mmol/L SR/ ms]

% Passive RyR leak - includes CaMKII regulation of leak flux
kleak = (1/3 + 10*RyR_CKp/3)*kleak;
J_SRleak = kleak*(y(31)-y(36));              %   [mmol/L cyt/ ms]

global Jleak_store 
Jleak_store(tStep,1) = J_SRCarel*Vsr/Vmyo + J_SRleak;   % Total Jleak [mmol/L cyt/ms]
Jleak_store(tStep,2) = J_SRleak;                        % Passive SR leak only [mmol/L cyt/ms]  
%% SERCA model - SR Ca uptake fluxes

% CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
fCKII_PLB = (1-.5*PLB_CKp);
%fracPKA_PLBo = .9926;       % Derived quantity - ((PLBtot - PLBp(baseline))/PLBtot)
fracPKA_PLBo = 1-0.007331;   % Derived quantity - ((PLBtot - PLBp(baseline))/PLBtot)
fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*3/4 + 1/4; % 0.25 with max PKA phosphorylation

% Select smaller value (resulting in max reduction of Kmf)
if fCKII_PLB < fPKA_PLB,
    Kmf = Kmf*fCKII_PLB;
else
    Kmf = Kmf*fPKA_PLB;
end

J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);

global Jserca_store  
Jserca_store(tStep) = J_serca;
%% Myofilament

% input contractile module ODE (6 state vars)
TSCa=y(54); TSCa_star=y(55); TSCa_tilde=y(56); 
TS_star=y(57); L_p=y(58); L_w=y(59); 

% PKA-dependent myofilament phosphoregulation (100 nM ISO)
fracPKA_Myoo = 0.003131; % Derived quantity (Myo_PKAp(baseline)/Myotot)
fracPKA_Myoiso = 0.9276; % Derived quantity (Myo_PKAp(ISO)/Myotot)
kPKA_Myo=(Myo_PKAp-fracPKA_Myoo)/(fracPKA_Myoiso-fracPKA_Myoo);
% Titin - Parallel elasticity
    uMyo=1; % no ISO effect w/ uMyo=0;
Ke=(1+uMyo*(0.5-1)*kPKA_Myo)*Ke; % 0.5* w/ ISO
% Crossbridges sensitivity (XBCa)
    uXBCa=1; % no ISO effect w/ uXBCa=0;
Zb=(1+uXBCa*(4.2-1)*kPKA_Myo)*Zb; % 4.2* w/ ISO
Zr=(1+uXBCa*(1.8-1)*kPKA_Myo)*Zr; % ...
Yr=(1+uXBCa*(2.2-1)*kPKA_Myo)*Yr;
% Crossbridges cycling (XBcy)
    uXBcy=1; % no ISO effect w/ uXBcy=0;
Za=(1+uXBcy*(1.24-1)*kPKA_Myo)*Za;
f=(1+uXBcy*(1.24-1)*kPKA_Myo)*f;
RLa=(1+uXBcy*(0.4-1)*kPKA_Myo)*RLa;
Zp=(1+uXBcy*(2.2-1)*kPKA_Myo)*Zp;
Yp=(1+uXBcy*(2.2-1)*kPKA_Myo)*Yp; 
Bp=(1+uXBcy*(3.4-1)*kPKA_Myo)*Bp;
Bw=(1+uXBcy*(3.4-1)*kPKA_Myo)*Bw;
Yc=(1+uXBcy*(0.4-1)*kPKA_Myo)*Yc;
Yd=(1+uXBcy*(2.2-1)*kPKA_Myo)*Yd;
Yv=(1+uXBcy*(1.6-1)*kPKA_Myo)*Yv;

% mechFlag defined for isometric (0) or isotonic (1) contraction
Liso=1.05; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if t>150 && t<=170,
%     Liso = Liso-(t-150)*(0.05/20);
% elseif t>170,
%     Liso = Liso-(170-150)*(0.05/20);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lm=1.05*mechFlag+Liso*(1-mechFlag);
Fm=0.87*mechFlag+0.87*(1-mechFlag);  
con=0;
corL=10;
L=1.03;
while abs(corL)>.00001,
   if mechFlag==0,
       Fm=alfa*(exp(bet*(Lm-L))-1);
   end
   FB=Ap*(TSCa_star+TS_star)*(L-L_p)+Aw*TSCa_tilde*(L-L_w); 
   w=FB+Ke*(L-Lz)^5+Le*(L-Lz)-Fm;
   w1=Ap*(TSCa_star+TS_star)+Aw*TSCa_tilde+5*Ke*(L-Lz)^4+Le+bet*(Fm+alfa);
   corL=-w/w1;
   L=L+.1*corL; 
   con=con+1;
   if con>200,
       break;
   end
end %while corL;   
TS=TSt-TSCa-TSCa_star-TSCa_tilde-TS_star;
ER=exp(-RLa*(L-La)^2); 
Yh=Yv*(1-exp(-gama*(L-L_w-hwr)^2));
if (L-L_w)>hwr;
    Yh=Fh*Yh;
end
fa=f*ER;
ga=Za+Yh;
gd=Yd*exp(-Yc*(L-Lc)); 
dTSCa=ga*TSCa_tilde-fa*TSCa+Yb*TS*y(38)^nc-Zb*TSCa;
dTSCa_star=Zr*TS_star*y(38)^nc-Yr*TSCa_star+Yp*TSCa_tilde-Zp*TSCa_star;
dTSCa_tilde=Zp*TSCa_star-Yp*TSCa_tilde+fa*TSCa-ga*TSCa_tilde;
dTS_star=-gd*TS_star+Yr*TSCa_star-Zr*TS_star*y(38)^nc;
dL_p=Bp*(L-L_p-hpr);
dL_w=Bw*(L-L_w-hwr);
if mechFlag==1,
    Lm=L+log((Fm+alfa)/alfa)/bet;
end 

% output contractile module ODE (dTSCa dTSCa_star dTSCa_tilde -> Ca buffer)
ydot(54:59)=[dTSCa dTSCa_star dTSCa_tilde dTS_star dL_p dL_w];

global Lmyo_store Fmyo_store  
Lmyo_store(tStep) = Lm;
Fmyo_store(tStep) = Fm;
%% Sodium and Calcium Buffering

% PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
fracTnIpo = .003131;  % Derived quantity (TnI_PKAp(baseline)/TnItot)
fPKA_TnI = (1.45-0.45*(1-TnI_PKAp)/(1-fracTnIpo)); % 1.45 with maximal phosphorylation
koff_tncl = koff_tncl*fPKA_TnI;

% PKA/CaMKII-dependent phosphoregulation of SERCA (NEW)
% Select smaller value (as for Kmf, see SERCA module)
if fCKII_PLB < fPKA_PLB,
    koff_sr = koff_sr*fCKII_PLB;
else
    koff_sr = koff_sr*fPKA_PLB;
end

% Na Buffers
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
%ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(19) = nc*(dTSCa+dTSCa_star+dTSCa_tilde); % myofilament
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = 0;% *** commented b/c buffering done by CaM module 
%ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
%J_CaB_cytosol = sum(ydot(19:25)); % CaB - error!
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);
%% Ion concentrations

% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);        % Csqn      [mM/ms]
ydot(31) = J_serca*Vmyo/Vsr-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30); % Ca_sr     [mM/ms] %Ratio 3 leak current

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   %[uA/uF]
ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
  +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec] 

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp;     % [uA/uF]
ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
  -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
  + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
%ydot(37)=0;
ydot(38) = -J_serca-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38)); % Cai

if Ca_clamp == 1,
    ydot(36) = 0; ydot(37) = 0; ydot(38) = 0; 
end

% junc_sl=J_ca_juncsl/Vsl*(y(36)-y(37));
% sl_junc=J_ca_juncsl/Vjunc*(y(37)-y(36));
% sl_myo=J_ca_slmyo/Vsl*(y(38)-y(37));
% myo_sl=J_ca_slmyo/Vmyo*(y(37)-y(38));

global Cai_store
Cai_store(tStep) = y(38);
%% Simulation type

% AP Waveform for AP clamp
global AP_Em AP_t

switch protocol
    case {'none',''}, % none
        I_app = 0;
        
    case 'pace', % pace w/ current injection at cycleLength 'cycleLength'
        if mod(t,cycleLength) <= 5
            I_app = 9.5;
        else
            I_app = 0.0;
        end

  	case 'v_step', % voltage step to 'input_par'
        step_period = 1e3;  % 1 Hz
        V_test = input_par;
        V_hold1 = -80; T_hold1 = 500;
        V_hold2 = V_test; T_hold2 = 4000;
	    if mod(t,step_period) <= T_hold1
            V_clamp = V_hold1;
        elseif mod(t,step_period) > T_hold1 && mod(t,step_period) <= T_hold1+T_hold2
            V_clamp = V_hold2;
        else
            V_clamp = V_hold1;
        end
		R_clamp = 0.01;
		I_app = (V_clamp-y(39))/R_clamp;
        
    case 'AP_clamp'
        % Determine appropriate voltage by interpolating between data points
        APc_period = 1e3; % 1 Hz
        ind1 = find(AP_t <= mod(t,APc_period),1,'last');
        ind2 = find(AP_t >= mod(t,APc_period),1,'first');
        tint = [AP_t(ind1),mod(t,APc_period),AP_t(ind2)];
        APint = interp1(AP_t,AP_Em,tint);
        potential = APint(2);
        
        R_clamp = .01;
        I_app = (potential-y(39))/R_clamp;   
               
end  
%% Membrane Potential

I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;                 % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk+Icftr;                         % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
ydot(39) = -(I_tot-I_app);

global Vmax_store
Vmax_store(tStep) = ydot(39);