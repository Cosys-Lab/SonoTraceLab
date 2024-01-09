function [mat_out] = normLog(mat_in,threshDB)
%NORMLIN Summary of this function goes here
%   Detailed explanation goes here
    thresh = 10^(threshDB/20);
    mat_out = mat_in / max(max(max(max(mat_in))));
    mat_out( mat_out < 0 ) = 0;
    mat_out = 20*log10( mat_out + thresh);
end