%% Full specific ideal-gas mixture entropy [J/(kg*K)]
% Input: temperature (T)[K], pressure (P)[Pa], mass fraction (Y), species SpS, molar masses (Mi) [kg/mol]
% Output: Mix entropy (S)[J/(kg*K)]
function S = mixStotal(T, P, Y, SpS, Mi)

global Runiv Pref

X = (Y ./ Mi) / sum(Y ./ Mi);
S = 0;

for i = 1:numel(SpS)
    if Y(i) > 0
        Ri = Runiv / Mi(i);
        S = S + Y(i) * ...
            (SNasa(T, SpS(i)) - Ri * log(X(i) * P / Pref));
    end
end

end