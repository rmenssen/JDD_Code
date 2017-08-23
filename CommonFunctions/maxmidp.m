function [maxmid, maxbin] = maxmidp(Ni, ri)

[~, maxbin] = max(Ni);
maxmid = ri(maxbin);

end