function RGB = gammaToneMapping(irrR, irrG, irrB, gamma)

Emax = max([irrR; irrG; irrB], [], 'all');
Emin = min([irrR; irrG; irrB], [], 'all');

Enorm(:,:,1) = (irrR - Emin) / (Emax - Emin);
Enorm(:,:,2) = (irrG - Emin) / (Emax - Emin);
Enorm(:,:,3) = (irrB - Emin) / (Emax - Emin);

RGB = Enorm .^ gamma;

end

