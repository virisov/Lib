function [sigW, sigI, tetS, epsW, epsI, tC] = loadMieSigma(freq,qq,MieDir)
%   function load precalculated Sigma matricies
%   and interpolate them to the required frequency 

if nargin<3
    MieDir = '';
end

sigW = [];
sigI = [];
tetS = [];
epsW = [];
epsI = [];

tC = -5;
T = 273+tC;

freqS = [5,(10:10:250)];
idx = find((freq>=freqS(1:end-1)) & (freq<freqS(2:end)));
if isempty(idx) 
    disp('Error: frequency is beyond limits 20-250GHz');
    return
end
f0 = freqS(idx);
f1 = freqS(idx+1);
x = log([f0,freq,f1]);
p = (x(3)-x(2))/(x(3)-x(1));

%medKey = 'w';   % w - water, i - ice

dName = sprintf('mie%03d',round(f0));

fName0 = [MieDir,dName,'/',dName,'_w',num2str(T)];
fName = sprintf('%s_q%02d.mat',fName0,round(10*qq));
if ~exist(fName,'file')
    disp(['Error: cannot find file ',fName,' in Mie library']);
    return
end
s0w = load(fName);
fName0 = [MieDir,dName,'/',dName,'_i',num2str(T)];
fName = sprintf('%s_q%02d.mat',fName0,round(10*qq));
if ~exist(fName,'file')
    disp(['Error: cannot find file ',fName,' in Mie library']);
    return
end
s0i = load(fName);

dName = sprintf('mie%03d',round(f1));

fName0 = [MieDir,dName,'/',dName,'_w',num2str(T)];
fName = sprintf('%s_q%02d.mat',fName0,round(10*qq));
if ~exist(fName,'file')
    disp(['Error: cannot find file ',fName,' in Mie library']);
    return
end
s1w = load(fName);
fName0 = [MieDir,dName,'/',dName,'_i',num2str(T)];
fName = sprintf('%s_q%02d.mat',fName0,round(10*qq));
if ~exist(fName,'file')
    disp(['Error: cannot find file ',fName,' in Mie library']);
    return
end
s1i = load(fName);

%nn = size(s0w.sigW);
%np =  prod(nn);
%sigW = reshape((reshape(s0w.sigW,np,1)*p + reshape(s1w.sigW,np,1)*(1-p),nn);
%sigI = reshape((reshape(s0i.sigW,np,1)*p + reshape(s1i.sigW,np,1)*(1-p),nn);

sigW = s0w.sigW*p + s1w.sigW*(1-p);
epsW = s0w.epsW*p + s1w.epsW*(1-p);

sigI = s0i.sigW*p + s1i.sigW*(1-p);
epsI = s0i.epsW*p + s1i.epsW*(1-p);

tetS = s0w.tet;

end



