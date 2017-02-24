function m = editParFile(vars,vals,fName0,fName1)
Nvars = length(vars);
f0 = fopen(fName0,'r');
f1 = fopen(fName1,'w');
tline = fgetl(f0);
m = 0;
while ischar(tline)
    for k=1:Nvars
        if ~isempty(strfind(tline,vars{k}))
            tline = [vars{k},' = ',vals{k},';     % new value'];
            m = m + 1;
        end
    end
    fprintf(f1,'%s\r\n',tline);
%    disp(tline);
    tline = fgetl(f0);
end
fclose(f0);
fclose(f1);
end