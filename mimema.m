function z = mimema(y,dim)
if nargin<2
    z = [min(y), mean(y), max(y)];
else
    z = [min(y,[],dim), mean(y,dim), max(y,[],dim)];
end
end

