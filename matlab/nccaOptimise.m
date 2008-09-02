function model = nccaOptimise(model,display,iters)

% NCCAOPTIMISE Optimise a given NCCA model
% FORMAT
% DESC model = nccaOptimise(model,display,iters) optimise a NCCA
% model
% ARG model : the model to be optimised
% ARG display : whether or not to display while optimisation
% proceeds, set to 2 for the most verbose and 0 for the least
% verbose.
% ARG iters : maximum number of iterations
% RETURN model : the optimised model
%
% SEEALSO : modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2007

% NCCA

if(nargin<3)
  iters = 1500;
  if(nargin<2)
    display = false;
    if(nargin<1)
      error('To Few Arguments');
    end
  end
end

if(isfield(model,'gy')&&~isempty(model.gy))
  model.gy = gpOptimise(model.gy,display,iters);
end
  
if(isfield(model,'fy')&&~isempty(model.fy))
  model.fy = gpOptimise(model.fy,display,iters);
end

if(isfield(model,'gz')&&~isempty(model.gz))
  model.gz = gpOptimise(model.gz,display,iters);
end
  
if(isfield(model,'fz')&&~isempty(model.fz))
  model.fz = gpOptimise(model.fz,display,iters);
end

% dynamics
tmp{3} = display;tmp{4} = iters;
if(isfield(model,'dyn_y')&&~isempty(model.dyn_y))
  model.dyn_y = modelOptimise(model.dyn_y,[],[],display,iters);
end

if(isfield(model,'dyn_z')&&~isempty(model.dyn_z))
  model.dyn_z = modelOptimise(model.dyn_z,[],[],display,iters);
end

return

