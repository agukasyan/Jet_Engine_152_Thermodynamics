%% Temperature dependent entropy J/(kg*K)
% Inputs: temperature (T)[K], mass fraction (Y), species SpS
% Output: temperature entropy (sT) [J/kg*k]
function sT = mixS(T, Y, SpS)

sT = zeros(size(T));

for i = 1:numel(SpS)
    sT = sT + Y(i) * SNasa(T, SpS(i));
end

end