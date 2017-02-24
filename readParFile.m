function par = readParFile(parFName)
% ------------------------------------------------
if ~exist(parFName,'file')
    par = [];
    return
end
fid = fopen(parFName);
tline = fgetl(fid);
while ischar(tline)
    disp(tline);
    eval(tline(1:end-1));
    tline = fgetl(fid);
end
fclose(fid); 
% ------------------------------------------------
end
