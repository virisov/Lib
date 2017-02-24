function j0 = getIndex(z,z0)
%   return the index on the nearest to z0 element of z
[~,j0] = min(abs(z - z0));
end

