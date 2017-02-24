function F = makeFreqStruct2(fSet)
%   fSet[NC,5] array:
%   1 - channel center frequency, GHz   
%   2 - number of side bands: 1 or 2   
%   3 - low freq.(2bands) or 0(1band), GHz   
%   4 - high freq.(2bands) or bandwidth(1band),  GHz
%   5 - antenna beam width(3dB), deg.
%
NC = size(fSet,1);
freq = [];
wf = [];
idxAnt = zeros(NC,1);
dAnt = [];
kdx = zeros(NC,2);
freqC = zeros(NC,1);
n = 1;
for i=1:NC
    ff = fSet(i,1);
    freqC(i) = ff;
    fLow = fSet(i,3);
    fHigh = fSet(i,4);
    fCent = (fLow+fHigh)/2;
    if fSet(i,2)==1
        m = 3;
        freq = [freq; ff-fHigh; ff; ff+fHigh];
        wf = [wf; [1;4;1]/6];
    else
        m = 6;
        freq = [freq; ff-fHigh; ff-fCent; ff-fLow; ff+fLow; ff+fCent; ff+fHigh];
        wf = [wf; [1;4;1;1;4;1]/12];
    end
    kdx(i,1) = n;
    kdx(i,2) = m;
    n = n+m;
    
    jdx = find(fSet(i,5)==dAnt);
    if isempty(jdx)
        dAnt = [dAnt, fSet(i,5)];
        idxAnt(i) = length(dAnt);
    else
        idxAnt(i) = jdx;
    end
        
end
F.NC = NC;      % number of channels
F.freqC = freqC;% center frequency of the channel
F.N = length(freq); % number of line frequencies: 3 or 6 per channel for 1 or 2-band
F.freq = freq;  % array of frequencies, GHz
F.wf = wf;      % weights for Simpson integration over bands 
F.kdx = kdx;    % kdx(1:NC,1) = n - start index in freq.array for a given channel, 
                % kdx(1:NC,2) = 3 or 6 number of freq. per channel
F.idxAnt = idxAnt;  % index of beam width: beam = dAnt(idxAnt(nc)), where nc - channel
F.dAnt = dAnt;  % antenna beam width array
