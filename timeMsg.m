%
% creates small string containing current time in HH:MM:SS format

function s = timeMsg ()

    c=uint8(clock);    
    s = sprintf('%02d:%02d:%02d', c(4), c(5), c(6));
    
end