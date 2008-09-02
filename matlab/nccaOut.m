function out = nccaOut(model,X,type,varargin)

% NCCAOUT Evaluate the output of an NCCA model.
% FORMAT
% DESC Y = nccaOut(model, x, type, ..) evaluates the output of a
% given NCCA model
% ARG model : the model for which the output is being evaluated.
% ARG x : the input position for which the output is required.
% ARG type : inference type
% 'YtoX' evaluate the latent cordinates from input Y
% 'ZtoX' evaluate the latent cordinates from input Z
% 'YtoZ' evaluate the corrsponding observations to Y
% 'ZtoY' evaluate the corrsponding observations to Z
% 'YtoXdyn' disambiguate the latent cordinates from  Y using
% dynamics
% 'ZtoXdyn' disambiguate the latent cordinates from Z using
% dynamics
% 'YtoZdyn' disambiguate the corresponding observations from Y
% using dynamics
% 'ZtoYdyn' disambiguate the corresponding observations from Z
% using dynamics
% ARG P1, P2,... : optional arguments to be passed
% RETURN : output of chose type
%
% SEEALSO : nccaCreate
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2007

% NCCA

if(nargin<3)
  if(model.fy.d==size(X,2))
    type = 'YtoZ';
  else
    type = 'ZtoY';
  end
  if(nargin<2)
    error('To Few Arguments');
  end
end

if(~strcmp(model.type,'ncca'))
  error('Wrong Model Type');
end

switch type
 case {'YtoZ','YtoX','YtoXdyn','YtoZdyn'}
  Y = X; clear X;
  if(length(varargin)<1)
    nr_nn = 10;
  end
  X = nccaComputeModes(model.gy,model.fz,Y,nr_nn,true);
  switch type
   case 'YtoX'
    out = X;
    return;
   case 'YtoZ'
    for(i = 1:1:size(X,1));
      out{i} = modelOut(model.fz,X{i});
    end
    return
   case {'YtoZdyn','YtoXdyn'}
    if(isfield(model,'dyn_z')&&~isempty(model.dyn_z))
      if(length(varargin)<3)
	iters = 100;
	if(length(varargin)<2)
	  balancing = 1;
	else
	  balancing = varargin{2};
	end
      else
	balancing = varargin{2};
	iters = varargin{3};
      end
      X = nccaComputeViterbiPath(model.gy,model.fz,model.dyn_z,X,Y,balancing,true);
      X = nccaSequenceOptimise(model.gy,model.fz,model.dyn_z,X,'independent',iters,balancing,true);
      switch type
       case 'YtoXdyn'
	out = X;
	return
       case 'YtoZdyn'
	out = modelOut(model.fz,X);
	return
      end
    else
      warning('No Dynamical Model Given: Returning Modes');
      out = X;
      return
    end
  end
 case {'ZtoY','ZtoX','ZtoXdyn','ZtoYdyn'}
  Z = X;clear X;
  if(length(varargin)<1)
    nr_nn = 10;
  end
  X = nccaComputeModes(model.gz,model.fy,Z,nr_nn,true);
  switch type
   case 'ZtoX'
    out = X;
    return
   case 'ZtoY'
    for(i = 1:1:size(X,1))
      out{i} = modelOut(model.fy,X{i});
    end
    return
   case {'ZtoXdyn','ZtoYdyn'}
    if(isfield(model,'dyn_y')&&~isempty(model.dyn_y))
      if(length(varargin)<3)
	iters = 100;
	if(length(varargin)<2)
	  balancing = 1;
	else
	  balancing = varargin{2};
	end
      else
	balancing = varargin{2};
	iters = varargin{3};
      end
      X = nccaComputeViterbiPath(model.gz,model.fy,model.dyn_y,X,Z,balancing,true);
      X = nccaSequenceOptimise(model.gz,model.fy,model.dyn_y,X,'independent',iters,balancing,true);
      switch type
       case 'ZtoXdyn'
	out = X;
	return
       case 'ZtoYdyn'
	out = modelOut(model.fy,X);
	return
      end
    else
      warning('No Dynamical Model Give: Returning Modes');
      out = X;
      return
    end
  end
 otherwise
  error('Unkown Type');
end




