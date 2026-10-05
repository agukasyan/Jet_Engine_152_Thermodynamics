%% Specific mixture gas constant [J/(kg*K)]
% Inputs: mass fraction (Y), molar masses (Mi) [kg/mol]
% Output: mixture gas constant (Rg) [J/kg*k]
function Rg = mixRg(Y, Mi)

global Runiv

Mmix = 1 / sum(Y ./ Mi);
Rg = Runiv / Mmix;

end