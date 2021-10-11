function RGB = reinhardGlobal(irrR, irrG, irrB, a, saturation)

% http://users.eecs.northwestern.edu/~ollie/eecs395/HW4/HW4.htm
% https://www.cl.cam.ac.uk/~rkm38/pdfs/mantiuk09cctm.pdf

Emax = max([irrR; irrG; irrB], [], 'all');
Emin = min([irrR; irrG; irrB], [], 'all');

Enorm(:,:,1) = (irrR - Emin) / (Emax - Emin);
Enorm(:,:,2) = (irrG - Emin) / (Emax - Emin);
Enorm(:,:,3) = (irrB - Emin) / (Emax - Emin);

L = rgb2gray(Enorm);

delta = 1e-6;
% Offset to avoid singularity
logLoffset = log(L + delta);
% Log average/geometric mean, formula in original paper seems to be
% slightly wrong
Lavg = exp(mean(logLoffset(:)));

T = (a / Lavg) .* L;
Tmax = max(T(:));

Ltone = T .* (1 + T ./ Tmax^2) ./ (1 + T);

% (2) in https://www.cl.cam.ac.uk/~rkm38/pdfs/mantiuk09cctm.pdf
for ch = 1:3
    RGB(:,:,ch) = (Enorm(:,:,ch) ./ L).^saturation .* Ltone;
end

end

