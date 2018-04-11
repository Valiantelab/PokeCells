function [lv] = locVar(T, locs, ap)
% Computes the Lv of Shinomoto

if ~isempty(locs)
    if length(locs) <= ap.io.lv_minisi
        lv = -1;
        return;
    end
    sum = 0.0;
    ISI = diff(T(locs));
    for i=1:length(ISI)-1
        sum = sum + ((ISI(i)-ISI(i+1))/(ISI(i)+ISI(i+1)))^2;
    end
    lv = 3/(length(ISI)-1)*sum;
else
    lv = -1;
end