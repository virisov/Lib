function [Garr,Gamma,gag] = makeGammas(zPTQW,GammaA,GammaS,F,A)
%
z = zPTQW(:,1); % altitude, km
P = zPTQW(:,2); % pressure, mb
T = zPTQW(:,3); % temperature, K
Q = zPTQW(:,4); % water vapor dencity, g/m^3
W = zPTQW(:,5); % liquid water content, g/m^3   
% W accounts for attenuation on small particles with negligible scattering

Nz = length(z);
Gamma0 = zeros(Nz,F.N);
for jf=1:F.N
   fr = repmat(F.freq(jf),Nz,1);
   Gamma0(:,jf) = abs_O2( T, P, Q, fr)  + abs_N2( T, P, fr);
   Gamma0(:,jf) = Gamma0(:,jf) + abs_H2O(T, P, Q, fr);
   Gamma0(:,jf) = Gamma0(:,jf) + abs_liquidH2O(T, W, fr);
   % O2 + N2 + WV absorption + small particles liquid (nep/km)
end
Gamma0 = Gamma0 + GammaA;   % + liquid/ice absorption
Gamma = Gamma0 + GammaS;    % + attenuation due to scattering  

%   ancilliary arrays needed for RTE integration in iterRTE2.m

G1 =  (Gamma(1:Nz-1,:) + Gamma(2:Nz,:))/2;
G = zeros(Nz-1,A.NDn,F.N);
dz = diff(z);
for jf=1:F.N
    G(:,:,jf) = (G1(:,jf).*dz)*A.secDn;
end
eG = exp(-G);
F1 = (G - 1 + eG)./G;
F2 = (1 - (1 + G).*eG)./G;
gag = Gamma0./Gamma;
gag1 =  (gag(1:Nz-1,:) + gag(2:Nz,:))/2;

GaG = zeros(Nz-1,A.NDn,F.N);
GaGHor = zeros(Nz,1,F.N);
for jf=1:F.N
    GaG(:,:,jf) = repmat(gag1(:,jf),1,A.NDn);
    GaGHor(:,1,jf) = gag(:,jf);
end

Garr.F1 = F1;
Garr.F2 = F2;
Garr.GaG = GaG;
Garr.GaGHor = GaGHor;
Garr.eG = eG;

end