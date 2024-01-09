function [mat_out] = normLin(mat_in)
%NORMLIN Summary of this function goes here
%   Detailed explanation goes here
    if(any(mat_in))
        mat_out = mat_in / max(max(max(max(mat_in))));
    else
        mat_out = mat_in;
    end
end