function [TupV, TupH] = surfaceScattSea(TdnV,TdnH,F,A,Tsurf)
%   function calculates reflected up-welling radiation
%   from flat sea water surface
%
%   Inputs:
%   TdnV, TdnH -- 2-D arrays of downwelling Tb at sea level (K)
%   1st index -- angle: A.tetDn (deg)
%   2nd index -- frequency: F.freq (GHz)
%   F - structure of frequencies, must have field F.freq (GHz)
%   A - structure of angles, must have A.tetDn field (deg)
%   Tsurf - surface temperature (K)
%
%   Outputs:
%   TupV, TupH -- 2-D arrays of upwelling reflected Tb at sea level (K)

sal = 35.0;
tempC = Tsurf - 273;

TupV = zeros(size(TdnV));
TupH = zeros(size(TdnH));
for jf = 1:F.N
    eps = ddE(sal,tempC,F.freq(jf));
    r = Fren2(eps,A.tetDn'*pi/180);
    TupV(:,jf) = Tsurf*(1-r(:,1)) + TdnV(:,jf).*r(:,1);
    TupH(:,jf) = Tsurf*(1-r(:,2)) + TdnH(:,jf).*r(:,2);
end

end
    