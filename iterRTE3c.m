function [Tv,Th,TauV,TauH] = iterRTE3c(Tv,Th,T,Beta,jdB,A,F,Garr,surfaceScatt,parSurf)
%
%   function runs an iteration of RTE
%
Tsky = 2.7;
Nz = size(Tv,1);
NDist = size(Beta,5);

%   particle scattering function

TauV = zeros(Nz,A.N,F.N);
TauH = TauV;
for jd=1:NDist
    iz = find(jdB==jd);
    if ~isempty(iz)
        for jf=1:F.N
            TauV(iz,:,jf) = ...
                Tv(iz,:,jf)*Beta(:,:,1,jf,jd)' + Th(iz,:,jf)*Beta(:,:,2,jf,jd)';
            TauH(iz,:,jf) = ...
                Tv(iz,:,jf)*Beta(:,:,3,jf,jd)' + Th(iz,:,jf)*Beta(:,:,4,jf,jd)';
        end
    end
end

%   downwelling profile calculation

TvDn = zeros(Nz,A.NDn,F.N);
TvDn(Nz,:,:) = Tsky;
ThDn = TvDn;

tv = TauV(:,(1:A.NDn),:);
th = TauH(:,(1:A.NDn),:);
for iz = Nz-1:-1:1
    p =  Garr.GaG(iz,:,:).*(T(iz).*Garr.F1(iz,:,:) + T(iz+1).*Garr.F2(iz,:,:));
    qv = (1 - Garr.GaG(iz,:,:)).*(tv(iz,:,:).*Garr.F1(iz,:,:) + tv(iz+1,:,:).*Garr.F2(iz,:,:));
    qh = (1 - Garr.GaG(iz,:,:)).*(th(iz,:,:).*Garr.F1(iz,:,:) + th(iz+1,:,:).*Garr.F2(iz,:,:));
    TvDn(iz,:,:) = TvDn(iz+1,:,:).*Garr.eG(iz,:,:) + p + qv;
    ThDn(iz,:,:) = ThDn(iz+1,:,:).*Garr.eG(iz,:,:) + p + qh;
end

%   surface scattering

[TvUp0,ThUp0] = surfaceScatt(squeeze(TvDn(1,:,:)),squeeze(ThDn(1,:,:)),F,A,parSurf);

%   upwelling profile calculation

TvUp = zeros(Nz,A.NDn,F.N);
ThUp = TvUp;
TvUp(1,:,:) = TvUp0;
ThUp(1,:,:) = ThUp0;

tv = TauV(:,(A.N:-1:A.NDn+2),:);
th = TauH(:,(A.N:-1:A.NDn+2),:);
for iz = 1:Nz-1
    p = Garr.GaG(iz,:,:).*(T(iz+1).*Garr.F1(iz,:,:) + T(iz).*Garr.F2(iz,:,:));
    qv = (1 - Garr.GaG(iz,:,:)).*(tv(iz+1,:,:).*Garr.F1(iz,:,:) + tv(iz,:,:).*Garr.F2(iz,:,:));
    qh = (1 - Garr.GaG(iz,:,:)).*(th(iz+1,:,:).*Garr.F1(iz,:,:) + th(iz,:,:).*Garr.F2(iz,:,:));
    TvUp(iz+1,:,:) = TvUp(iz,:,:).*Garr.eG(iz,:,:) + p + qv;
    ThUp(iz+1,:,:) = ThUp(iz,:,:).*Garr.eG(iz,:,:) + p + qh;
end

%   horizontal look

ja = A.NDn+1;
TvHor = Garr.GaGHor.*repmat(T,[1,1,F.N]) + (1 - Garr.GaGHor).*TauV(:,ja,:);
ThHor = Garr.GaGHor.*repmat(T,[1,1,F.N]) + (1 - Garr.GaGHor).*TauH(:,ja,:);

%   merge complete angle scan

Tv = cat(2,TvDn,TvHor);
Tv = cat(2,Tv,flipdim(TvUp,2));
Th = cat(2,ThDn,ThHor);
Th = cat(2,Th,flipdim(ThUp,2));

end