function [P,T] = standardAtmosphere(z)
% standard atmosphere
% z, m
% P, Pa
% T, K

G = 9.80655;    % m/s^2
R = 287.04;     % m^2/(K*s^2)

T0 = 288.15;    % 15C
P0 = 101325.0;  % Pa

h11 = 11000;    % m
P11 = 22632.0;  % Pa
T11 = 216.65;

j11 = z>h11;
T = T0 - 0.0065*z;
T(j11) = T11;
P = P0*(1 - 0.0065*z/T0).^5.2561;
P(j11) = P11*exp(-G*(z(j11)-h11)./(R*T11));

end

