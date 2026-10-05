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

