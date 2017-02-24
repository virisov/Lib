function A = angleGrid(DTet)
%
%  uniform angle grid (0 - zenith)
%  over-zenith and over-nadir angle range added for
%  subsequent antenna convolution
%
%D = 400*eps;

A.DTet = DTet;   % deg.
A.tetDn = (0:DTet:90-DTet);  % deg. 0 - zenith
A.secDn = 1./cos(A.tetDn*pi/180);
A.tetD = [A.tetDn, 90, 180 - fliplr(A.tetDn)];
A.NDn = length(A.tetDn);
A.N = length(A.tetD);
end
