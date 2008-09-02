function options = nccaOptions(approx)

% NCCAOPTIONS Return default options for NCCA model.
% FORMAT
% DESC options = nccaOptions(approx) returns the default options in
% a structure for a NCCA model.
% ARG approx : approximation type, either 'ftc' (no approximation),
% 'dtc' (deterministic training conditional), 'fitc' (fully
% independent training conditional) or 'pitc' (partially
% independent training conditional).
% RETURN options : option structure
%
% SEEALSO : nccaCreate, gpOptions
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2007

% NCCA

if(nargin<1)
  approx = 'ftc';
end

options = gpOptions(approx);
options.optimiser = 'scg';
options.scale2var1 = true;
options.kern = {'cmpnd','rbfard','linard','bias','white'};

return