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
% Ideal adiabatic diffuser; negligible outlet velocity

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

    % Air enthalpy at the trial temperature
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

%% Compressor [2-3]
% Ideal adiabatic compressor; negligible inlet and outlet velocities

v3 = 0;

% Outlet pressure [Pa].
P3 = P3overP2 * P2;

% Target temperature entrophy [J/(kg*K)]
s3target = mixS(T2, Yair, SpS) + RgAir * log(P3/P2);

TL = T2;                    % Lower temperature bound [K]
TH = 3000;                  % Upper temperature bound [K]

% Sanity check
if mixS(TL, Yair, SpS) > s3target || ...
        mixS(TH, Yair, SpS) < s3target
    error('Compressor temperature is outside the search interval.');
end

while (TH - TL) > 0.01
    Ti = (TL + TH)/2;

    % Temperature part of air entropy at the trial temperature
    s3trial = mixS(Ti, Yair, SpS);

    if s3trial > s3target
        TH = Ti;            % Trial temperature is too high
    else
        TL = Ti;            % Trial temperature is too low
    end
end

T3 = (TL + TH)/2;

h3 = mixH(T3, Yair, SpS);    % Outlet enthalpy [J/kg]
wc = h3 - h2;               % Specific work input [J/kg air]
Wc = mair * wc;             % Compressor power input [W]

% Display сompressor results
fprintf('Compressor: T3 = %.2f K, P3 = %.2f kPa\n', ...
    T3, P3/1000);
fprintf('Compressor power = %.3f MW\n', Wc/1e6);


%% Combustor [3-4]: composition
% Complete combustion with excess oxygen
% Species order: fuel, O2, CO2, H2O, N2

% Fuel element composition [O H C N Ar]
fuelElements = SpS(1).Elcomp;
zO = fuelElements(1);
yH = fuelElements(2);
xC = fuelElements(3);

% [mol O2/mol fuel]
nuO2 = xC + yH/4 - zO/2;

% Reaction coefficients
nu = [-1, -nuO2, xC, yH/2, 0];

% Combined air and fuel mass fractions before combustion
Y3 = (mair * Yair + mfurate * Yfuel) / mtot;

% Convert mass flow rates to molar flow rates [mol/s]
n3 = mtot * Y3 ./ Mi;

% Complete combustion of the incoming fuel [mol/s]
n4 = n3 + nu * n3(1);

n4(n4 < 0) = 0;              % Remove tiny negative round-off values

% Product mass fractions
Y4 = n4 .* Mi / sum(n4 .* Mi);

% Mass conservation residual [kg/s]
massError = sum(n4 .* Mi) - mtot;

% Gas constants before and after combustion [J/(kg*K)]
Rg3 = mixRg(Y3, Mi);
Rg4 = mixRg(Y4, Mi);

% Stoichiometric air-fuel ratio and equivalence ratio
AFstoich = nuO2 * MAir / (Xair(2) * Mi(1));
phi = AFstoich / AF;


%% Combustor [3-4]: energy balance
% Adiabatic, no shaft work, negligible kinetic energy
% No pressure loss; fuel inlet temperature = T3

P4 = P3;
v4 = 0;
Tfuel = T3;

% Incoming enthalpy flow [W]
hfuel = mixH(Tfuel, Yfuel, SpS);
H3in = mair * h3 + mfurate * hfuel;

% Required specific enthalpy of products [J/kg]
% NASA enthalpies include formation enthalpy
h4 = H3in / mtot;

% Find T4 using bisection
TL = T3;
TH = 3000;

if mixH(TL, Y4, SpS) > h4 || ...
        mixH(TH, Y4, SpS) < h4
    error('Combustor temperature is outside the search interval.');
end

while (TH - TL) > 0.01
    Ti = (TL + TH)/2;
    h4trial = mixH(Ti, Y4, SpS);

    if h4trial > h4
        TH = Ti;
    else
        TL = Ti;
    end
end

T4 = (TL + TH)/2;

% Display сombustor results
fprintf('Combustor: T4 = %.2f K, P4 = %.2f kPa, phi = %.4f\n', ...
    T4, P4/1000, phi);


%% Turbine [4-5]
% Ideal adiabatic turbine driving the compressor.

v5 = 0;

% Shaft power balance [J/kg]
h5 = h4 - Wc / mtot;

% Find outlet temperature using bisection
TL = 200;
TH = T4;

if mixH(TL, Y4, SpS) > h5 || ...
        mixH(TH, Y4, SpS) < h5
    error('Turbine temperature is outside the search interval.');
end

while (TH - TL) > 0.01
    Ti = (TL + TH)/2;

    if mixH(Ti, Y4, SpS) > h5
        TH = Ti;
    else
        TL = Ti;
    end
end

T5 = (TL + TH)/2;

% Outlet pressure from the isentropic condition [Pa]
P5 = P4 * exp((mixS(T5, Y4, SpS) ...
    - mixS(T4, Y4, SpS)) / Rg4);

fprintf('Turbine: T5 = %.2f K, P5 = %.2f kPa\n', ...
    T5, P5/1000);


%% Nozzle [5-6]
% Ideal adiabatic nozzle expanding to ambient pressure

P6 = Pamb;

if P5 <= P6
    error('Nozzle inlet pressure must exceed ambient pressure.');
end

% Target temperature part of entropy
s6target = mixS(T5, Y4, SpS) + Rg4 * log(P6/P5);

% Find outlet temperature using bisection
TL = 200;
TH = T5;

if mixS(TL, Y4, SpS) > s6target || ...
        mixS(TH, Y4, SpS) < s6target
    error('Nozzle temperature is outside the search interval.');
end

while (TH - TL) > 0.01
    Ti = (TL + TH)/2;

    if mixS(Ti, Y4, SpS) > s6target
        TH = Ti;
    else
        TL = Ti;
    end
end

T6 = (TL + TH)/2;
h6 = mixH(T6, Y4, SpS);

% Outlet velocity from energy conservation [m/s]
v6 = sqrt(2 * (h5 - h6) + v5^2);

fprintf('Nozzle: T6 = %.2f K, P6 = %.2f kPa, v6 = %.2f m/s\n', ...
    T6, P6/1000, v6);
%% Results: state arrays and Table 1
T = [T1 T2 T3 T4 T5 T6];
P = [P1 P2 P3 P4 P5 P6];
v = [v1 v2 v3 v4 v5 v6];

% States 1-3 contain air; states 4-6 contain combustion products.
Ycomp = {Yair, Yair, Yair, Y4, Y4, Y4};

h = zeros(1, 6);
S = zeros(1, 6);

for k = 1:6
    h(k) = mixH(T(k), Ycomp{k}, SpS);
    S(k) = mixStotal(T(k), P(k), Ycomp{k}, SpS, Mi);
end

fprintf('\nTable 1: thermodynamic states\n');
fprintf('%-14s', 'State');
fprintf('%11d', 1:6);
fprintf('\n');

fprintf('%-14s', 'P [kPa]');
fprintf('%11.2f', P/1000);
fprintf('\n');

fprintf('%-14s', 'T [K]');
fprintf('%11.2f', T);
fprintf('\n');

fprintf('%-14s', 'v [m/s]');
fprintf('%11.2f', v);
fprintf('\n');

fprintf('%-14s', 'h [kJ/kg]');
fprintf('%11.2f', h/1000);
fprintf('\n');

fprintf('%-14s', 's [kJ/(kg K)]');
fprintf('%11.4f', S/1000);
fprintf('\n');

%% Results: Table 2
fprintf('\nTable 2: combustor mass fractions\n');
fprintf('AF = %.2f, AFstoich = %.4f, phi = %.4f\n', ...
    AF, AFstoich, phi);

fprintf('%-14s %12s %12s\n', 'Species', 'Initial', 'Final');

for i = 1:numel(SpS)
    fprintf('%-14s %12.6f %12.6f\n', ...
        cSpecies{i}, Y3(i), Y4(i));
end

fprintf('%-14s %12.3f %12.3f\n', ...
    'Rg [J/(kg K)]', Rg3, Rg4);
