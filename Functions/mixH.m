%% Specific mixture enthalpy J/kg 
% Inputs: temperature (T) [K], mass fractions (Y), species (SpS)
% Output: Enthalpy (h) [J/kg] 
function h = mixH(T, Y, SpS)

h = zeros(size(T));

for i = 1:numel(SpS)
    h = h + Y(i) * HNasa(T, SpS(i));
end

end