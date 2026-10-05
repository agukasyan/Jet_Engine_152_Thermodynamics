clear all;close all;clc;
addpath('General');
addpath('Functions');


%% 152 conditions
v1 = 200;          % m/s
Tamb = 250;        % K
P3overP2 = 7;
Pamb = 55000;      % Pa
mfurate = 0.68;    % kg/s
AF = 102.78;       % kg air / kg fuel
cFuel = 'CH4';

%% Nasa implementation
load(fullfile('General', 'NasaThermalDatabase.mat'));

% Constants NASA
global Runiv Pref
Runiv = 8.314472;    % J/(mol*K)
Pref = 1.01235e5;   % Pa

% NASA database
cSpecies = {cFuel,'O2','CO2','H2O','N2'};
iSp = myfind({Sp.Name}, cSpecies);
SpS = Sp(iSp);
Mi = [SpS.Mass];


%% Air calculations
% Air mole fractions: 21% O2 and 79% N2.
Xair = [0 0.21 0 0 0.79];
MAir = Xair * Mi';
Yair = Xair .* Mi / MAir;
Yfuel = [1 0 0 0 0];

% Specific gas constant of air J/(kg*K)
RgAir = mixRg(Yair, Mi);

% Mass flow rates
mair = AF * mfurate;
mtot = mair + mfurate;


%% Diffuser [1-2]
% Ideal adiabatic diffuser; negligible outlet velocity.

T1 = Tamb;
P1 = Pamb;
v2 = 0;

NSp = length(SpS);
Rg = RgAir;                  % Air gas constant [J/(kg*K)]

% Inlet enthalpy
for i = 1:NSp
    hi(i) = HNasa(T1, SpS(i));
end

h1 = Yair * hi';             % Air enthalpy at state 1 [J/kg]

% Energy conservation
h2 = h1 + 0.5*v1^2 - 0.5*v2^2;

% Find T2 using bisection
TL = T1;                    % Lower temperature bound [K]
TH = 1000;                  % Upper temperature bound [K]
iter = 0;

while abs(TH - TL) > 0.01
    iter = iter + 1;
    Ti = (TL + TH)/2;

    % Air enthalpy at the trial temperature.
    for i = 1:NSp
        hi2(i) = HNasa(Ti, SpS(i));
    end

    h2i = Yair * hi2';

    if h2i > h2
        TH = Ti;            % Trial temperature is too high
    else
        TL = Ti;            % Trial temperature is too low
    end
end

T2 = (TH + TL)/2;

% Properties at the inlet and outlet
for i = 1:NSp
    hi2(i) = HNasa(T2, SpS(i));
    si1(i) = SNasa(T1, SpS(i));
    si2(i) = SNasa(T2, SpS(i));
end

h2check = Yair * hi2';       % Evaluated enthalpy at the found T2
s1thermal = Yair * si1';
s2thermal = Yair * si2';

% Outlet pressure from the isentropic condition
lnPr = (s2thermal - s1thermal)/Rg;
P2 = P1 * exp(lnPr);

% Check the entropy change at fixed air composition
S1 = s1thermal - Rg*log(P1/Pref);
S2 = s2thermal - Rg*log(P2/Pref);

% Display diffuser results
fprintf('Diffuser: T2 = %.2f K, P2 = %.2f kPa\n', T2, P2/1000);
fprintf('Enthalpy residual = %.6f J/kg\n', h2check - h2);