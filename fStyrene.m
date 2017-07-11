function dydL=etano(L,y);
nj=y(1:10); % [mol/s] molar flowrates of each component
T=y(11);    % [K] temperature 
P=y(12);    % [Pa] pressure
global vis A  R E ko Href alfa dQdL Tref Cpj  
global MW dp Por d0  dcat  

Xj=nj/sum(nj)    ; % molar   fraction of each component, [1x10]
Pj=Xj*P          ; % [Pa]    partial pressure of each  component, [1x10]
PM=MW*Xj         ; % [g/mol] average molecular weight
Qv=sum(nj)*R*T/P ; % [m3/s]  volumetric flow

d=(P/R/T)*PM / 10^3;  % [kg/m3] density of the gas mixture, [1x1] 

H=Href +Cpj*alfa'*(T-Tref); % [J/mol] enthalpy of reaction to the temperature, [1x6]
k=ko.*exp(-E/(R*T))       ; % [mol/m3/s] constant rate, [1x6]
K = exp(-1.6) * 101325    ; % [Pa] equilibrium constant of the 1st reaction, [1x1]

%  Reactions kinetics
%  E-benz  Styrene   H2  Ethylene  Toluene   CO    CO2     CH4    H2O     Benz
r(1)=k(1)* ( Pj(1) - Pj(2)*Pj(3)/K) ;
r(2)=k(2)*Pj(1)                     ;
r(3)=k(3)*Pj(1)*Pj(3)               ;
r(4)=k(4)*Pj(4)*Pj(9)^2             ; 
r(5)=k(5)*Pj(8)*Pj(9)               ;
r(6)=k(6)*Pj(6)*Pj(9)               ;

dcat = 2350 ; % [kg/m3] catalyst density, [1x1]
u = Qv/A    ; % [m/s] gas velocity though reactor, [1x1]
G = d * u   ; % [kg/m2/s] mass flowrate divided by area, [1x1]


% System of differential equations:

dydL(1:10) =  A * r * alfa                                           ;  %  Molar balance

dydL(11)   = (dQdL-A*(r*H'))/(Cpj*nj)                                ;  % Energy balance


dydL(12)   = - 100 * ( G * (1-Por) / (d0  *dp * Por^3 ) ) * ...
                ( 150 * (1 - Por) * vis / dp + 1.75 * G )            ;  % Momentum Balance (Ergun Eq.)
            
            
dydL=dydL' ;



