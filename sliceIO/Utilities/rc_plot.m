function [r c] = rc_plot(nplots)

if nplots < 4
    r = nplots;
    c = 1;
    return
end

nsq = sqrt(nplots);

if (floor(nsq))^2 == nplots
    % Its a perfect square
    r = nsq;
    c = nsq;
    return;
end

r = round(nsq);
c = r;

if (r*c) < nplots
    c = c + 1;
end