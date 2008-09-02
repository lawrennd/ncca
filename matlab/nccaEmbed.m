function [Xsy Xsz Xy Xz] = nccaEmbed(Y,Z,data_var_keep,shared_var_keep,consolidation_var_keep,verbose)

% NCCAEMBED Consolidate two data sets using NCCA
% FORMAT
% DESC  [Xsy Xsz Xy Xz] =
% nccaEmbed(Y,Z,data_var_keep,shared_var_keep,consolidation_var_keep,VERBOSE)
% creates a NCCA model structure with default parameter settings as
% specified by the options vector.
% ARG Y : Observations in space 1 (either in form on Kernel or
% Datapoints (in which case linear algorithm will be applied)
% ARG Z : Observations in space 2
% ARG data_var_keep : Percentage of variance to keep after
% PCA-reduction
% ARG shared_var_keep : Minimum percentage of cannonical variate
% shared dimension needs to explain.
% ARG consolidation_var_keep : Percentage of variance of data
% Consolidation needs to explain
% ARG VERBOSE
% RETURN Xsy : Shared latent representation in Y space
% RETURN Xsz : Shared latent representation in Z space
% RETURN Xy : Independent latent representation of Y
% RETURN Xz : Independent latent representation of Z
%
% SEEALSO : NCCACREATE, NCCAOPTIMISE
%
% COPYRIGHT : Neil D. Lawrence, Carl Henrik Ek, 2007

% NCCA

if(nargin<6)
  verbose = false;
  if(nargin<5)
    consolidation_var_keep = 10;
    if(nargin<4)
      shared_var_keep = 50;
      if(nargin<3)
	data_var_keep = 99;
	if(nargin<2)
	  error('To Few Arguments');
	end
      end
    end
  end
end

% check input arguments
if(size(Y,1)~=size(Z,1))
  error('Different number of observations in Y and Z');
end

if(size(Y,1)==size(Y,2))
  if(verbose)
    fprintf(';-----------------------------------;\n');
    fprintf('Y given in terms of Kernel Matrix\n');
  end
  Ky = Y;
  clear Y;
else
  if(verbose)
    fprintf(';-----------------------------------;\n');
    fprintf('Y given in terms of points\n');
  end
  Ky = cov(Y');
  clear Y;
end
if(size(Z,1)==size(Z,2))
  if(verbose)
    fprintf('Z given in terms of Kernel Matrix\n');
    fprintf(';-----------------------------------;\n');
  end
  Kz = Z;
  clear Z;
else
  if(verbose)
    fprintf('Z given in terms of Points\n');
    fprintf(';-----------------------------------;\n');
  end
  Kz = cov(Z');
  clear Z;
end

% pre-process Kernels
Ky = kernelCenter(Ky);
Kz = kernelCenter(Kz);
Ky = Ky./sum(diag(Ky));
Kz = Kz./sum(diag(Kz));
Ky = (Ky+Ky')./2;
Kz = (Kz+Kz')./2;

% Re-represent data in terms of dominant Principal Directions
[lambda_y,Vy] = eigdec_valid(Ky);
[lambda_z,Vz] = eigdec_valid(Kz);

lambda_y = lambda_y./sum(lambda_y);
lambda_z = lambda_z./sum(lambda_z);
[lambda_y,ind] = sort(lambda_y,'descend');
Vy = Vy(:,ind);
[lambda_z,ind] = sort(lambda_z,'descend');
Vz = Vz(:,ind);

cut_off = data_var_keep/100;
c_lambda_y = cumsum(lambda_y);
c_lambda_z = cumsum(lambda_z);

nr_feature_basis_y = length(find(c_lambda_y<cut_off))+1;
nr_feature_basis_z = length(find(c_lambda_z<cut_off))+1;
clear c_lambda_y;
clear c_lambda_z;
clear cut_off;

if(verbose)
  fprintf(';-----------------------------------;\n');
  fprintf('PCA reduction:\n');
  fprintf('Y:\t %dD\n',nr_feature_basis_y);
  fprintf('Z:\t %dD\n',nr_feature_basis_z);
  fprintf(';-----------------------------------;\n');
end

% KPCA Expansion Coefficients
Ay = Vy*inv(sqrtm(diag(lambda_y)));
Ay = Ay(:,1:1:nr_feature_basis_y);
Az = Vz*inv(sqrtm(diag(lambda_z)));
Az = Az(:,1:1:nr_feature_basis_z);

% Reduced Feature space representation
Cy = Ky*Ay;
Cz = Kz*Az;

% design matrix in reduced feature space
Cyy = Cy'*Cy;
Czz = Cz'*Cz;

% Set-Up CCA
Cyy = (Cyy+Cyy')./2;
Czz = (Czz+Czz')./2;

Cyy_inv = pdinv(Cyy);
Czz_inv = pdinv(Czz);

By = Cyy_inv*Cy'*Cz*Czz_inv*Cz'*Cy;
Bz = Czz_inv*Cz'*Cy*Cyy_inv*Cy'*Cz;
% CCA
[lambda_y,alpha_y] = eigdec_valid(By);
[lambda_z,alpha_z] = eigdec_valid(Bz);

lambda_y = lambda_y./sum(lambda_y);
lambda_z = lambda_z./sum(lambda_z);

[lambda_y,ind] = sort(lambda_y,'descend');
alpha_y = alpha_y(:,ind);
[lambda_z,ind] = sort(lambda_z,'descend');
alpha_z = alpha_z(:,ind);


% Re-represent shared space with dominant cannonical components
cut_off = shared_var_keep/100;
nr_shared_basis = length(find(lambda_y>cut_off));
if(nr_shared_basis==0)
  nr_shared_basis = 1;
  if(verbose)
    fprintf(';-----------------------------------;\n');
    fprintf('No Shared Basis Explains %3.0f of the shared variance\n',shared_var_keep);
    fprintf('Forcing the Shared Dimensionality to one in order to proceed\n');
    fprintf(';-----------------------------------;\n');
  end
end
alpha_y = alpha_y(:,1:1:nr_shared_basis);
alpha_z = alpha_z(:,1:1:nr_shared_basis);

% NCCA in PCA space
beta_y = [];
beta_z = [];
for(i = 1:1:nr_feature_basis_y-nr_shared_basis)
  beta_y(:,i) = ncca_lin(Cyy,[alpha_y beta_y],true);
end
for(i = 1:1:nr_feature_basis_z-nr_shared_basis)
  beta_z(:,i) = ncca_lin(Czz,[alpha_z beta_z],true);
end

% enforce orthogonality
%tmp = cgrscho([alpha_y beta_y]);
%alpha_y = tmp(:,1:1:size(alpha_y,2));
%beta_y = tmp(:,size(alpha_y,2)+1:1:end);
%tmp = cgrscho([alpha_z beta_z]);
%alpha_z = tmp(:,1:1:size(alpha_z,2));
%beta_z = tmp(:,size(alpha_z,2)+1:1:end);
%clear tmp;

% Directions in Kernel induced or Feature Space
%alpha_y = Ay*alpha_y;
%alpha_z = Az*alpha_z;
%beta_y = Ay*beta_y;
%beta_z = Az*beta_z;

% enforce orthogonality
%tmp = cgrscho([alpha_y beta_y]);
%alpha_y = tmp(:,1:1:size(alpha_y,2));
%beta_y = tmp(:,size(alpha_y,2)+1:1:end);
%tmp = cgrscho([alpha_z beta_z]);
%alpha_z = tmp(:,1:1:size(alpha_z,2));
%beta_z = tmp(:,size(alpha_z,2)+1:1:end);
%clear tmp;

% normalise
alpha_y = normalise_vector(alpha_y);
alpha_z = normalise_vector(alpha_z);
beta_y = normalise_vector(beta_y);
beta_z = normalise_vector(beta_z);

% reduce dimensionality of independent spaces
Ny_shared = sum(diag(alpha_y'*Cyy*alpha_y))/trace(Cyy);
Nz_shared = sum(diag(alpha_z'*Czz*alpha_z))/trace(Czz);

lambda_y = diag(beta_y'*Cyy*beta_y);
lambda_z = diag(beta_z'*Czz*beta_z);

c_lambda_y = cumsum(lambda_y)+Ny_shared;
c_lambda_z = cumsum(lambda_z)+Nz_shared;

cut_off = (consolidation_var_keep/(data_var_keep/100))/100;
nr_basis_independent_y = min(length(find(c_lambda_y<cut_off))+1,length(c_lambda_y));
nr_basis_independent_z = min(length(find(c_lambda_z<cut_off))+1,length(c_lambda_z));

beta_y = beta_y(:,1:1:nr_basis_independent_y);
beta_z = beta_z(:,1:1:nr_basis_independent_z);

Ny_independent = sum(diag(beta_y'*Cyy*beta_y))/trace(Cyy);
Nz_independent = sum(diag(beta_z'*Czz*beta_z))/trace(Czz);

if(verbose)
  fprintf(';-----------------------------------;\n');
  fprintf('Consolidation:\n');
  fprintf('Shared dimension\t\t %3.0f\n',nr_shared_basis);
  fprintf('Y independent dimension\t\t %3.0f\n',nr_basis_independent_y);
  fprintf('Z independent dimension\t\t %3.0f\n',nr_basis_independent_z);
  fprintf(';-----------------------------------;\n');
  fprintf('Variance Modeling:\n');
  fprintf('Shared:\n');
  fprintf('\t\tY: %2.1f\tZ: %2.1f\n',Ny_shared*100,Nz_shared*100);
  fprintf('Independent:\n');
  fprintf('\t\tY: %2.1f\tZ: %2.1f\n',Ny_independent*100,Nz_independent*100);
  fprintf(';-----------------------------------;\n');
end

% Embeddings
Xsy = Ky*(Ay*alpha_y);
Xsz = Kz*(Az*alpha_z);
Xy = Ky*(Ay*beta_y);
Xz = Kz*(Az*beta_z);

return

%--------------------------------------------------------%
function beta = ncca_lin(K,alpha,verbose)

if(nargin<3)
  verbose = false;
end

A = K - alpha*alpha'*K;
[D beta] = eigdec_valid(A,1);

if(verbose)
  if(any(~isreal(D)))
    warning('NCCA: Complex Eigenvalues');
  end
  if(length(find(D<0))~=0)
    warning('NCCA: Negative Eigenvalues');
  end
end

return

function v = normalise_vector(v)

for(i = 1:1:size(v,2))
  v(:,i) = v(:,i)./norm(v(:,i),2);
end

return


function [lambda V] = eigdec_valid(A,dim,zero_trunc)

if(nargin<3)
  zero_trunc = 0;
  if(nargin<2)
    dim = [];
    if(nargin<1)
      error('To Few Arguments');
    end
  end
end

[lambda,V] = eigdec(A,size(A,1));

% remove imaginary solutions
if(~isreal(lambda))
  ind = find(~imag(lambda));
  lambda = lambda(ind);
  V = V(:,ind);
end

% normalise and sort
lambda = lambda./sum(lambda);
[lambda ind] = sort(lambda,'descend');
V = V(:,ind);

% remove "zero" solutions
ind = find(abs(lambda)>zero_trunc);
lambda = lambda(ind);
V = V(:,ind);

% re-normalise and re-sort(not needed)
lambda = lambda./sum(lambda);
[lambda ind] = sort(lambda,'descend');
V = V(:,ind);

if(~isempty(dim))
  if(length(lambda)<dim)
    error('Problem does not contain %d valid dimensions\n',dim);
  end
  V = V(:,1:1:dim);
  lambda = lambda(1:1:dim);
end

return

