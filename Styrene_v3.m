clear all;clc;close all;

global vis A  R E ko Href alfa dQdL Tref Cpj  
global MW dp Por d0  dcat 

%% PROBLEM DATA
M=184.7 * 10^3 / 3600 ;  % [mol/s] total input molar flowrate feed

%  E-benz  Styrene       H2  Ethylene  Toluene   CO    CO2     CH4    H2O     Benz
MW=[   106     104.15    2    25.05    92.14    28.01  44.01  16.04   18      78.11   ];             % [g/mol] molecular weight
xi=[  0.0772   0.0001    0     0      0.0003      0     0       0    0.9224     0     ];             %  input molar fraction of feed
ni=[  14.25   1.34e-2   0     0      6.19e-2     0     0       0    170.42   0.0002  ] * 1000/3600 ; % [mol/s] input molar flowrate to reactor
      
alfa=[  -1        1      1     0          0       0     0       0      0        0   ;   %1  
        -1        0      0     1          0       0     0       0      0        1   ;   %2
        -1        0     -1     0          1       0     0       1      0        0   ;   %3
         0        0      4    -1          0       2     0       0     -2        0   ;   %4
         0        0      3     0          0       1     0      -1     -1        0   ;   %5
         0        0      1     0          0      -1     1       0     -1        0  ] ;  %6   stoichometric coefficients of reactions
     
Cpj=[301.795803874600 277.188875848083 30.2388208892725 91.9827389727751 250.975686103559 32.7070798588799 52.7026167287852 69.0339845402977 40.3392353348438 202.636986558307]; % [J/mol] molar heat capacity
    
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
Pen=1.6*10^5 ; % [Pa] input pressure to reactor
Ten=700+273  ; % [K]  input temperature to reactor
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                                                                                                                                 
                                                      
%     1                                                2                       3                                    4                              5                                 6                                                         
Href= [117.6                                         105.5                  -54.68                           210.11819921875                     206.124                         -41.386  ] * 10^3        ; % [J/mol] enthalpy of reaction
E =   [90792.7                                     207944.6                  91462.1                           103972.3                          65688.73                      73554.64   ]               ; % [J/mol] activation energy
ko=   [1967 * 10^3 / 101325                7.3e8 * 10^3 / 101325      1748 * 10^3 / (101325^2)          1209 * 10^3 / (101325^3)      69.11 * 10^3 / (101325^2)        7329 * 10^3 / (101325^2)     ]     ; % [mol/m3/Pa/s] pre-exponential factor

Radio=1.95/2 ; % [m] reactor radius
D=1.95       ; % [m] reactor diameter

%dQdL=80*2*pi()*Radio*10^3; %J/m s
dQdL=0  ;  % adiabatic process implies

Tref=298       ; % [K] reference temperature
R=8.314        ; % [J/mol/K] ideal gases constant
A = pi*Radio^2 ; % [m2] reactor area
a=0.5;
vis=2.86*10^-5 ; % [Pa·s] viscosity
dp= 1*10^-3    ; % [m] particle diameter
Por= 0.445     ; % porosity
PMi=MW.*xi     ; % [g/mol] average molecular weight, [1x1]
dcat=2350      ; % [kg/m3] catalys density

Qv0=sum(ni)*R*Ten/Pen ; % [m3/s]
d0 = sum(ni.*xi)/Qv0  ; % [kg/m3]
u0 = Qv0/A            ; % [m/s]       


%% solution to the problem
L = linspace(0,3,100); % [m] independient variable fixed
CI=[ni Ten Pen]      ; % [mol/s, K, Pa] initial conditions
lengthspan = L * Por ; % [m] length where there is reaction 

[Lcat Y]=ode23s(@fStyrene, lengthspan, CI) ;

njout=Y(:,1:10);        % [mol/s] molar flowrates along reactor
T=Y(:,11);              % [K]     temperature along reactor
P=Y(:,12);              % [Pa]    pressure alog reactor
njout=njout*3600/1000 ; % [kmol/h]

%FIG. MOLAR FLOWRATES OF ALL COMPONENTS AND STYRENE
figure(1)
subplot(1,2,1)
plot(L,njout)
xlabel('L (m)')
ylabel('nj (kmol/h)')
legend('E-benz','Styrene','H_2','Ethylene','Toluene','CO','CO2','CH4','H2O','Benz')
grid on
subplot(1,2,2)
plot(L,njout(:,2))
xlabel('L (m)')
ylabel('n styrene (kmol/h)')
grid on

%FIG. TEMPERATURE
T=T-273; %ºC
figure(2)
plot(L,T)
xlabel('L (m)')
ylabel('T (ºC)')
grid on
hold on
T1=[680.941088287670;667.245758583485;656.847200532591;648.615372963946;641.907873511842;636.327141314535;631.617429711722;627.569453757844;624.049104169784;620.955713662092;618.215379195658;615.769517329006;613.589340992974;611.618667152107;609.828129233792;608.193716856662;606.695622211755;605.317382913046;604.045237027503;602.867632747233];
L1=[0.0750000000000000;0.225000000000000;0.375000000000000;0.525000000000000;0.675000000000000;0.825000000000000;0.975000000000000;1.12500000000000;1.27500000000000;1.42500000000000;1.57500000000000;1.72500000000000;1.87500000000000;2.02500000000000;2.17500000000000;2.32500000000000;2.47500000000000;2.62500000000000;2.77500000000000;2.92500000000000];
plot(L1,T1,'r')
legend('T Matlab (ºC)','T Hysys (ºC)')

%plot(L1,T1)
%FIG. PRESSURE
figure(3)
plot(L,P/10^5)
xlabel('L (m)')
ylabel('P (bar)')
grid on
hold on
P1=[157.359067998848;154.741912820437;152.130227756791;149.512039504346;146.878847326925;144.224003962127;141.541909409230;138.827523239360;136.076285178511;133.283796688961;130.445687895271;127.557485503994;124.614520079665;121.611716463727;118.543720135908;115.404668841940;112.188057642302;108.886612387034;105.492082669494;101.995018246982];
L1=[0.0750000000000000;0.225000000000000;0.375000000000000;0.525000000000000;0.675000000000000;0.825000000000000;0.975000000000000;1.12500000000000;1.27500000000000;1.42500000000000;1.57500000000000;1.72500000000000;1.87500000000000;2.02500000000000;2.17500000000000;2.32500000000000;2.47500000000000;2.62500000000000;2.77500000000000;2.92500000000000];
plot(L1,P1/100,'r')
legend('T Matlab (ºC)','T Hysys (ºC)')


% FIG.  ¡E-BENCENE CONVERSION
conv = (ni(1)-njout(:,1)*1000/3600)/ni(1);
figure(4)
plot(L,conv)
xlabel('L (m)') 
ylabel('conv E-Benzene')
grid on

% FIG. CONCENTRACION
for i=1:100
    Qv_1(i) = sum(njout(i,:))*R*T(i)./P(i) ;  % [ m3/h ] caudal volumétrico a lo largo del reactor
    C_1(i,:) = njout(i,:)/Qv_1(i)          ;  % [ kmol/m3 ] concentraciones de los componentes a lo largo del reactor
end

figure(5)
plot(L,C_1)
xlabel('L (m)')
ylabel('C (mol/m3')
grid on
legend('E-benz','  Styrene','   H_2','    Ethylene','  Toluene ','  CO ','     CO2 ','     CH4 ','    H2O  ','       Benz')



% Extents of reactions calculations

extentReaction(1) =   njout(end,2) - njout(1,2)             ; %[ kmol/h ]
extentReaction(2) =  njout(end,10) - njout(1,10)            ; %[ kmol/h ]
extentReaction(3) =   njout(end,5) - njout(1,5)             ; %[ kmol/h ]

extentReaction(4) = - ( ( njout(end,4) - njout(1,4)   )  - ...
                             extentReaction(2)    )         ; %[ kmol/h ] 
                         
extentReaction(5) = - ( ( njout(end,8) - njout(1,8)   )  - ...
                             extentReaction(3)    )         ; %[ kmol/h ]
                         
extentReaction(6) = ( njout(end,7) - njout(1,7)   )         ; %[ kmol/h ]     



figure(6)
nEB=[14.2499999866417; 12.4560162424313;11.1711503037447;10.1967313614372;9.42819219494128;8.80432980811443;8.28655079725042;7.84909924249188;7.47421020013982;7.14909814158854;6.86430286263879;6.61264822868618;6.38858718208035;6.18770973533945;6.00655445577194;5.84232362108173;5.69273208734722;5.55589673547117;5.43025433311547;5.31449974481609;5.20753896091912];
nST=[0; 1.56641672472810;2.69608190077576;3.56004836013360;4.24469288154870;4.80189187780074;5.26493799450298;5.65632779029577;5.99173212923244;6.28250946809277;6.53711057617288;6.76197297406702;6.96208968773460;7.14143696940739;7.30314656441206;7.44975370996405;7.58333152148602;7.70558965675207;7.81794783709231;7.92159128920470;8.01751296275359];
L1=[0; 0.0750000000000000;0.225000000000000;0.375000000000000;0.525000000000000;0.675000000000000;0.825000000000000;0.975000000000000;1.12500000000000;1.27500000000000;1.42500000000000;1.57500000000000;1.72500000000000;1.87500000000000;2.02500000000000;2.17500000000000;2.32500000000000;2.47500000000000;2.62500000000000;2.77500000000000;2.92500000000000];
plot(L1,nEB,'r')
grid on
hold on
plot(L1,nST,'r')
xlabel('L (m)')
ylabel('nj (kmol/h)')

plot(L,njout(:,1),'b')
grid on
hold on
plot(L,njout(:,2),'b')
xlabel('L (m)')
ylabel('nj (kmol/h)')

legend('E-Benzene and Styrene flowrates in Matlab','E-Benzene and Styrene flowrates in Hysys')

