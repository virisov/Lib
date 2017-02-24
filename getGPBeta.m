function [beta,GammaA,GammaS] = getGPBeta(freq,A,Dist,MieDir)
%   function reads frequency interpolated (in loadIceSigma) 
%   Sig-scattering matricies, calculated by a chain of scripts
%   runGreenF.m, runSigma77.m, and rotateParticle77.m, 
%   and integrate them according to G-distribution of sizes
%   with parameters in Dist structure array
%   Size of Dist-array corresponds to number of different uniform layers
%   Dist.C - fraction of liquid particles in water/ice mixture
%   GammaA, GammaS - absorption and attenuation due to scattering, nep/km
%
%   n(a)*da = A*(Lambda*a)^(alpha-1)*exp(-Lambda*a)*da
%
%   a - particle radius, cm
%   Lambda, 1/cm
%   A, 1/cm^4 
%
if nargin<4
    MieDir = '';
end
wk = 2*pi*freq/30.0;	%	rad/cm
q = (0.1:0.1:10);  % Mie library particle size
if ~isfield(Dist,'P')   % just in case
    Dist.P = 0;
end
    
Ntet = A.N;
tet = A.tetD*pi/180;
stet2 = (2*pi*mean(diff(tet)))*repmat(sin(tet),Ntet,1);

ALk2 = Dist.ALam/(Dist.Lambda*wk^2);	% cm^-1	coeff
kLq6 = (wk/(Dist.Lambda*q(1)))^6;
kLq3 = (wk/(Dist.Lambda*q(1)))^3;
alfa = Dist.alpha;

betaW = zeros(Ntet,Ntet,4);
betaI = zeros(Ntet,Ntet,4);
betaPI = zeros(Ntet,Ntet,4);
gAPI = 0;   % create them to avoid error if Dist.P==0 and 
gSPI = 0;   % if-Dist.P block is skipped

dq = mean(diff(q));
dx = dq*Dist.Lambda/wk;
Lq = length(q);
for iq=1:Lq
    x = Dist.Lambda*q(iq)/wk;
    z = ALk2*x^(alfa - 1)*exp(-x)*dx;		% cm^-1

    [sigW, sigI, tetS, epsW, epsI] = loadMieSigma(freq,q(iq),MieDir);
    if isempty(sigW)    % cannot find file in Mie library
        beta = [];
        GammaA = [];
        GammaS = [];
        return
    end

    if iq==1

        %	Rayleigh approximation
        
        bs = intGamma0(alfa+6,x);   % numerical intGamma doesn't
        ba = intGamma0(alfa+3,x);   % work well at alfa=1;
        
        sw = sigW*ALk2*kLq6*bs;	% cm^-1
        si = sigI*ALk2*kLq6*bs;	% cm^-1
        [ScW,AbW] = smallSphereSc(epsW,q(1));   % used only for Gamma
        gSW = ScW*ALk2*kLq6*bs;
        gAW = AbW*ALk2*kLq3*ba;
        [ScI,AbI] = smallSphereSc(epsI,q(1));
        gSI = ScI*ALk2*kLq6*bs;
        gAI = AbI*ALk2*kLq3*ba;
        
        %   trapec.sum.rule
        
        sw = sw + sigW*z/2;
        si = si + sigI*z/2;
        gSW = gSW + ScW*z/2;
        gAW = gAW + AbW*z/2;
        gSI = gSI + ScI*z/2;
        gAI = gAI + AbI*z/2;
        
    else

        %	Mie solution
        
        sw = sw + sigW*z;
        si = si + sigI*z;
        [S,ScW,AbW] = SphereSc2(sqrt(epsW),q(iq),0);    % used only for Gamma
        [S,ScI,AbI] = SphereSc2(sqrt(epsI),q(iq),0);
        gSW = gSW + ScW*z;
        gAW = gAW + AbW*z;
        gSI = gSI + ScI*z;
        gAI = gAI + AbI*z;
    end
end

sw = sw*1e5;  % 1/km
si = si*1e5;  % 1/km
%   in general case of non-spherical particles gSW-gAI and
%   GammaA, GammaS depend on polarization
gSW = gSW*1e5;	%	1/km
gAW = gAW*1e5;	%	1/km
gSI = gSI*1e5;	%	1/km
gAI = gAI*1e5;	%	1/km

[tx1,ty1] = meshgrid(tet);
dtetS = mean(diff(tetS));
%   make grid to go through the poles for 
%   smooth interpolation over them
[tx0,ty0] = meshgrid([tetS(1)-dtetS, tetS, tetS(end)+dtetS]);

for ip=1:4
    z = [sw(:,1,ip), sw(:,:,ip), sw(:,end,ip)];
    z = [z(1,:); z; z(end,:)];
    b = interp2(tx0,ty0,z,tx1,ty1,'spline');
    betaW(:,:,ip) = b.*stet2;
    z = [si(:,1,ip), si(:,:,ip), si(:,end,ip)];
    z = [z(1,:); z; z(end,:)];
    b = interp2(tx0,ty0,z,tx1,ty1,'spline');
    betaI(:,:,ip) = b.*stet2;
end

%   now add non-spherical particles
if Dist.P>eps
    qp = (0.1:0.1:3.0);
    LqP = length(qp);
    dq = mean(diff(q));     % the same
    dx = dq*Dist.Lambda/wk; % the same
    for iq=1:LqP
        x = Dist.Lambda*qp(iq)/wk;
        z = ALk2*x^(alfa - 1)*exp(-x)*dx;		% cm^-1

        [SigPI, AbPI, tetPS] = loadIceSigma(freq,qp(iq));
        
        if isempty(sigW)    % cannot find file in Particle library
            beta = [];
            GammaA = [];
            GammaS = [];
            return
        end
        if iq==1   %	Rayleigh approximation
            bs = intGamma0(alfa+6,x);   % numerical intGamma doesn't
            ba = intGamma0(alfa+3,x);   % work well at alfa=1;
            si = SigPI*ALk2*kLq6*bs;	% cm^-1
            gAPI = AbPI*ALk2*kLq3*ba;

            si = si + SigPI*z/2;
            gAPI = gAPI + AbPI*z/2;
        else
            si = si + SigPI*z;
            gAPI = gAPI + AbPI*z;
        end
    end
    si = si*1e5;  % 1/km
    gAPI = gAPI*1e5;	%	1/km    
    %   Achtung! 
    %   In general case of non-spherical particles gSW-gAI and
    %   GammaA, GammaS depend on polarization
    %   Needs to be done in future! (Add total absorption through P-vector)
    dtetS = mean(diff(tetPS));
    %   make grid to go through the poles for 
    %   smooth interpolation over them
    [tx0,ty0] = meshgrid([tetPS(1)-dtetS, tetPS', tetPS(end)+dtetS]);
    for ip=1:4
        z = [si(:,1,ip), si(:,:,ip), si(:,end,ip)];
        z = [z(1,:); z; z(end,:)];
        b = interp2(tx0,ty0,z,tx1,ty1,'spline');
        betaPI(:,:,ip) = b.*stet2;
    end
    bs = squeeze(sum(betaPI,2));
    gSPI2 = [mean(bs(:,1) + bs(:,2)), mean(bs(:,3) + bs(:,4))]; % 2-pol
    gSPI = mean(gSPI2);   % for now until absorption polarization won't be resolved
end

%   Dist.P is a fraction of non-spherical ice particles among all ice
%   fraction (1 - Dist.C)

beta = betaW*Dist.C + (betaI*(1 - Dist.P) + betaPI*Dist.P)*(1 - Dist.C);
GammaA = gAW*Dist.C + (gAI*(1 - Dist.P) + gAPI*Dist.P)*(1 - Dist.C);
GammaS = gSW*Dist.C + (gSI*(1 - Dist.P) + gSPI*Dist.P)*(1 - Dist.C);

%   Ensure exact energy conservation
%   Test shows that energy deffect is less than 1%

bs = squeeze(sum(beta,2));
beta(:,:,1:2) = beta(:,:,1:2).*(GammaS/mean(bs(:,1) + bs(:,2)));
beta(:,:,3:4) = beta(:,:,3:4).*(GammaS/mean(bs(:,3) + bs(:,4)));

%   now we need to normalize 'beta' on 'GammaS' 
%   [see iterRTE2.m and RTE_algorithm.pdf. We get here 'beta_with_tilda' Eg.8-11]

beta = beta/GammaS;     % dimensionless

end

