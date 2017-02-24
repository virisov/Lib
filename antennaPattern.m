function [tt,dd] = antennaPattern(w05, tetMax, DTet)
%{
tetD = (0:0.25:10);
ka = [30 40 50 60 70 80];
col = 'bgrcmk';
figure, hold
for i=1:6
   plot(tetD,F1(ka(i),tetD),['.-',col(i)])
   tz(i) = fzero(@(x) F1(ka(i),x)-0.5 ,[0 10]);
end
grid

figure
for i=1:6
   semilogy(tetD,F1(ka(i),tetD),['.-',col(i)])
   if i==1
      hold
   end
end
grid

p = polyfit(tz,1./ka,2);
tt = [1:0.1:4];
v = polyval(p,tt);

figure, hold
plot(tz,ka,'o');
plot(tt,1./v,'-');
%}

p = [-0.00000346525727   0.00840738763534  -0.00000701877567]; % 3<w05<8
v = polyval(p,w05/2);
ka = 1/v;
t = (0:DTet:tetMax);
d = F1(ka,t);
tt = [-fliplr(t(2:end)), t];
dd = [fliplr(d(2:end)), d];
end

function f = F1(ka,tetD)
Z1 = 2.4048;
x = ka*sin(tetD*pi/180);
y = Z1^2*besselj(0,x)./(Z1^2 - x.^2);
f = y.^2;
end
