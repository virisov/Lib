function [TupV, TupH] = surfaceScattSeaLand(TdnV,TdnH,F,A,parSurf)
%   function calculates reflected up-welling radiation
%   from flat land-sea surface
%
%   Inputs:
%   TdnV, TdnH -- 2-D arrays of downwelling Tb at sea level (K)
%   1st index -- angle: A.tetDn (deg)
%   2nd index -- frequency: F.freq (GHz)
%   F - structure of frequencies, must have field F.freq (GHz)
%   A - structure of angles, must have A.tetDn field (deg)
%   parSurf = [Tsurf(K), landSeaMask(0-sea), salinity(0/00)]
%
%   Outputs:
%   TupV, TupH -- 2-D arrays of upwelling reflected Tb at sea level (K)

Tsurf = parSurf.Tsurf;
rWL = parSurf.rWL;

TupV = zeros(size(TdnV));
TupH = zeros(size(TdnH));
for jf = 1:F.N
    TupV(:,jf) = Tsurf*(1-rWL(:,1,jf)) + TdnV(:,jf).*rWL(:,1,jf);
    TupH(:,jf) = Tsurf*(1-rWL(:,2,jf)) + TdnH(:,jf).*rWL(:,2,jf);
end

end
    