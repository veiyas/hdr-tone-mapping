function RGB = gammaLuminanceToneMapping(irrR, irrG, irrB, gamma, sat)

Emax = max([irrR; irrG; irrB], [], 'all');
Emin = min([irrR; irrG; irrB], [], 'all');

Enorm(:,:,1) = (irrR - Emin) / (Emax - Emin);
Enorm(:,:,2) = (irrG - Emin) / (Emax - Emin);
Enorm(:,:,3) = (irrB - Emin) / (Emax - Emin);

L = rgb2gray(Enorm);
Ltone = L .^ gamma;

% (2) in https://www.cl.cam.ac.uk/~rkm38/pdfs/mantiuk09cctm.pdf
for ch = 1:3
    RGB(:,:,ch) = (Enorm(:,:,ch) ./ L).^sat .* Ltone;
end

end

