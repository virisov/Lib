function xF = averBand(x,F)
%
%   averaging over frequency band
%
nn = size(x);
xF = zeros(nn(1),nn(2),F.NC);
for jf = 1:F.NC
    for kf = 1:F.kdx(jf,2)
        m = F.kdx(jf,1) + kf - 1;
        xF(:,:,jf) = xF(:,:,jf) + x(:,:,m)*F.wf(m);
    end
end

end
