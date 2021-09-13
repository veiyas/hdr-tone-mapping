function value = weightingFunction(z)
% See (4) in "Recovering High Dynamic Range Radiance Maps from Photographs"
Zmin = 0;
Zmax = 255;
middle = (Zmin + Zmax) / 2;

if z <= middle
    value = z - Zmin;
else
    value = Zmax - z;
end

end

