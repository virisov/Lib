function [TupV, TupH] = surfaceScattBB(TdnV,TdnH,F,A,Tsurf)
%   function calculates reflected up-welling radiation
%   from black body surface
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

TupV = zeros(size(TdnV));
TupH = zeros(size(TdnH));
for jf = 1:F.N
    TupV(:,jf) = Tsurf;
    TupH(:,jf) = Tsurf;
end

end
    