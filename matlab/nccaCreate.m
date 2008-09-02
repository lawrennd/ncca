function model = nccaCreate(void1,void2,Y,Z,data_var_keep,shared_var_keep,consolidation_var_keep,shared_space,options)

% NCCACREATE Create a NCCA model
% FORMAT
% DESC 	MODEL = NCCACREATE(void, void,
% Y,Z,DATA_VAR_KEEP,SHARED_VAR_KEEP,CONSOLIDATION_VAR_KEEP,SHARED_SPACE,OPTIONS)
% creates a NCCA model structure with default parameter settings as
% specified by the options vector.
% ARG void : for compatibility
% ARG void : for compatibility
% ARG Y : Observations in space 1 (either in form on Kernel or
% Datapoints (in which case linear algorithm will be applied)
% ARG Z : Observations in space 2
% ARG data_var_keep : Percentage of variance to keep after
% PCA-reduction
% ARG shared_var_keep : Minimum percentage of cannonical variate
% shared dimension needs to explain
% ARG consolidation_var_keep : Percentage of variance of data
% Consolidation needs to explain
% ARG shared_space : Shared Space from 'Y' or 'Z'
% ARG options : options structure as defined by nccaOptions.m
% RETURN model : model structure containing the Gaussian process.
%
% SEEALSO : nccaOptions, nccaCreate, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2007

% NCCA

if(nargin<9)
  options = nccaOptions('ftc');
  if(nargin<8)
    shared_space = 'Z';
    if(nargin<7)
      consolidation_var_keep = 90;
      if(nargin<6)
	shared_var_keep = 10;
	if(nargin<5)
	  data_var_keep = 97;
	  if(nargin<4)
	    error('To Few Arguments');
	  end
	end
      end
    end
  end
end

clear void1 void2;

% learn Embedding
[Xsy Xsz Xy Xz] = nccaEmbed(Y,Z,data_var_keep,shared_var_keep,consolidation_var_keep,true);

switch shared_space
 case 'Y'
  Xs = Xsy;
  % align Z to Y
  [void void transf] = procrustes(Xsy,Xsz);
  %Xz = transf.b*Xz*transf.T + transf.c;
 case 'Z'
  Xs = Xsz;
  % align Y to Z
  [void void transf] = procrustes(Xsz,Xsy);
  %Xy = transf.b*Xy*transf.T + transf.c;
 otherwise
  Xs = Xsz;
  % align Z to Y
  [void void transf] = procrustes(Xsz,Xsy);
  %Xy = transf.b*Xy*transf.T + transf.c;
end

% create model
model.type = 'ncca';
model.N = size(Y,1);
model.ds = size(Xs,2);
model.dy = size(Xy,2);
model.dz = size(Xz,2);

model.gy = gpCreate(size(Y,2),size(Xs,2),Y,Xs,options);

if(~isempty(Xy))
  model.fy = gpCreate(size(Xs,2)+size(Xy,2),size(Y,2),[Xs Xy],Y,options);
else
  model.fy = [];
end

model.gz = gpCreate(size(Z,2),size(Xs,2),Z,Xs,options);

if(~isempty(Xz))
  model.fz = gpCreate(size(Xs,2)+size(Xz,2),size(Z,2),[Xs Xz],Z,options);
else
  model.fz = [];
end

return