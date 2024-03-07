function R=deg2rad(D)
%DEG2RAD Convert angles from degrees to radians
%
%  rad = DEG2RAD(deg) converts angles from degrees to radians.
%
%  See also RAD2DEG, DEG2DMS, ANGLEDIM, ANGL2STR

% Copyright 1996-2006 The MathWorks, Inc.
% Written by:  E. Byrns, E. Brown
% $Revision: 1.9.4.2 $    $Date: 2006/05/24 03:33:47 $


if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(D)
     warning('Imaginary parts of complex ANGLE argument ignored')
     D = real(D);
end

R = D*pi/180;
