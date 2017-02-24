function epsL = epsLand(N)
% N random realizations of eps
mse2 = [3.1895,    1.3360];
epsL = [];
while length(epsL)<N
    ra = randn(2*N,2);
    epsSa = mse2(1) + mse2(2)*ra(:,2);
    te = abs(0.1 + 0.02*ra(:,1));
    te(epsSa<1) = [];
    epsSa(epsSa<1) = [];
    epsL = [epsL; epsSa.*complex(1, te)./sqrt(1 + te.^2)];
end
epsL(N+1:end) = [];
end